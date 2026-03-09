#!/bin/bash
#PBS -P <PROJECT>
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l mem=1000GB
#PBS -l jobfs=500GB
#PBS -l ncpus=48
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

# ============================================================================
# Filter Plant Contamination and Extract Clean Microbial Reads
# ============================================================================
# Description: Uses competitive mapping results to separate plant-derived
#              reads from microbial reads, then creates clean FASTQ files
# Input: Sorted BAM files from competitive mapping, original FASTQ files
# Output: Clean microbial FASTQ files (SMAG + unmapped reads)
# ============================================================================

# Exit immediately on error and provide detailed logging
set -euo pipefail
set -x

# NOTE: Use ~half the allocated CPUs for each samtools step in pipes to avoid oversubscription.
THREADS=$(( ${PBS_NCPUS:-48} / 2 ))
[ "$THREADS" -lt 1 ] && THREADS=1

# Load necessary modules
source <CONDA_PATH>/etc/profile.d/conda.sh
conda activate <ALIGN_ENV_PATH>

module load samtools

# Change to script directory (where classify_reads.py is located)
cd "<SCRIPT_DIR>"

# --- Define variables ---
BAM_DIR="<INPUT_BAM_DIR>"
REF_LIST_DIR="<REFERENCE_CONTIG_LISTS_DIR>"
FASTQ_DIR="<ORIGINAL_FASTQ_DIR>"
FINAL_OUTPUT_DIR="${BAM_DIR}/final_output"

# Create necessary output directories
mkdir -p "$FINAL_OUTPUT_DIR"
mkdir -p "${BAM_DIR}/split_reads"
mkdir -p "${BAM_DIR}/idxstats"

# Main summary file for all barcode proportions
PROPORTIONS_SUMMARY_FILE="${FINAL_OUTPUT_DIR}/all_barcodes_read_proportions.txt"

# Overwrite or create the summary file and write the header.
# We'll dynamically generate the header to include all reference names.
HEADER="Barcode\tTotal_Mapped_Reads"
for ref_list in "$REF_LIST_DIR"/*_contigs.txt; do
    ref_name=$(basename "$ref_list" _contigs.txt)
    HEADER+="\t${ref_name}_reads(%)"
done
echo -e "$HEADER" > "$PROPORTIONS_SUMMARY_FILE"

# --- Main Loop to process each barcode ---
for INPUT_BAM in "$BAM_DIR"/*.sorted.bam; do
    
    # Get the barcode name from the filename
    BARCODE_NAME=$(basename "$INPUT_BAM" .sorted.bam)
    echo "--- Processing $BARCODE_NAME ---"

    # --- Task 1: Split reads using the Python script for advanced tie-breaking ---
    echo "  - Running Python script for read splitting..."
    SPLIT_OUTPUT_DIR="${BAM_DIR}/split_reads/${BARCODE_NAME}"
    mkdir -p "$SPLIT_OUTPUT_DIR"
    
    # Analyzes all alignments for each read (include primary and secondary alignments)
    python classify_reads.py "$INPUT_BAM" "$REF_LIST_DIR" "$SPLIT_OUTPUT_DIR"
    
    echo "  - Building exact BEDs from BAM header (idxstats -> BED)..."
    IDXSTATS_FILE="${BAM_DIR}/idxstats/${BARCODE_NAME}.idxstats.txt"
    samtools idxstats "$INPUT_BAM" > "$IDXSTATS_FILE"

    # AWK is safer than grep -w for contig names with dots/pipes/etc.
    TEMP_PLANT_BED=$(mktemp)
    TEMP_SMAG_BED=$(mktemp)

    awk 'NR==FNR{want[$1]=1; next} ($1 in want){print $1"\t0\t"$2}' \
        "${REF_LIST_DIR}/plant_contigs.txt" "$IDXSTATS_FILE" > "$TEMP_PLANT_BED"

    awk 'NR==FNR{want[$1]=1; next} ($1 in want){print $1"\t0\t"$2}' \
        "${REF_LIST_DIR}/smag_contigs.txt"  "$IDXSTATS_FILE" > "$TEMP_SMAG_BED"

    # Based on the read name lists, create the final BAM files
    # Add a second filtering step to ensure that only alignments to the designated
    # contigs are kept.
    echo "  - Creating final BAM files..."

    # NOTE: Do the filtering in two passes to keep memory usage flat.
    # 1) name filter -> temp BAM; 2) region filter -> final BAM. Then index.

    # PLANT
    if [ -s "${SPLIT_OUTPUT_DIR}/plant_unique_names.txt" ] && [ -s "$TEMP_PLANT_BED" ]; then
        samtools view -@ "$THREADS" -b -N "${SPLIT_OUTPUT_DIR}/plant_unique_names.txt" "$INPUT_BAM" \
          -o "${SPLIT_OUTPUT_DIR}/plant_unique_temp.bam"
        samtools index -c "${SPLIT_OUTPUT_DIR}/plant_unique_temp.bam"
        samtools view -@ "$THREADS" -b -L "$TEMP_PLANT_BED" \
          -o "${SPLIT_OUTPUT_DIR}/plant_unique_final.bam" \
          "${SPLIT_OUTPUT_DIR}/plant_unique_temp.bam"
        rm -f "${SPLIT_OUTPUT_DIR}/plant_unique_temp.bam" "${SPLIT_OUTPUT_DIR}/plant_unique_temp.bam.csi"
    else
        echo "    (plant) no names or no contigs; emitting empty BAM"
        samtools view -H "$INPUT_BAM" | samtools view -b -o "${SPLIT_OUTPUT_DIR}/plant_unique_final.bam"
    fi

    # SMAG
    if [ -s "${SPLIT_OUTPUT_DIR}/smag_unique_names.txt" ] && [ -s "$TEMP_SMAG_BED" ]; then
        samtools view -@ "$THREADS" -b -N "${SPLIT_OUTPUT_DIR}/smag_unique_names.txt" "$INPUT_BAM" \
          -o "${SPLIT_OUTPUT_DIR}/smag_unique_temp.bam"
        samtools index -c "${SPLIT_OUTPUT_DIR}/smag_unique_temp.bam"
        samtools view -@ "$THREADS" -b -L "$TEMP_SMAG_BED" \
          -o "${SPLIT_OUTPUT_DIR}/smag_unique_final.bam" \
          "${SPLIT_OUTPUT_DIR}/smag_unique_temp.bam"
        rm -f "${SPLIT_OUTPUT_DIR}/smag_unique_temp.bam" "${SPLIT_OUTPUT_DIR}/smag_unique_temp.bam.csi"
    else
        echo "    (smag) no names or no contigs; emitting empty BAM"
        samtools view -H "$INPUT_BAM" | samtools view -b -o "${SPLIT_OUTPUT_DIR}/smag_unique_final.bam"
    fi
    
    # AMBIGUOUS
    if [ -s "${SPLIT_OUTPUT_DIR}/ambiguous_names.txt" ]; then
        samtools view -@ "$THREADS" -N "${SPLIT_OUTPUT_DIR}/ambiguous_names.txt" -b "$INPUT_BAM" \
          -o "${SPLIT_OUTPUT_DIR}/ambiguous_final.bam"
    else
        echo "    (ambiguous) no names; emitting empty BAM"
        samtools view -H "$INPUT_BAM" | samtools view -b -o "${SPLIT_OUTPUT_DIR}/ambiguous_final.bam"
    fi

    echo "  - Indexing final BAM files..."
    # The -c flag is essential for large genomes to create a .csi index
    samtools index -c "${SPLIT_OUTPUT_DIR}/plant_unique_final.bam"
    samtools index -c "${SPLIT_OUTPUT_DIR}/smag_unique_final.bam"
    samtools index -c "${SPLIT_OUTPUT_DIR}/ambiguous_final.bam"
    echo "  - Task 1 finished for $BARCODE_NAME."

    # --- Task 2: Calculate proportions for each crop and SMAG ---
    echo "  - Starting Task 2: Calculating read proportions..."
    
    total_mapped_reads=$(awk '{sum+=$3} END {print sum+0}' "$IDXSTATS_FILE")
    
    proportions_line="$BARCODE_NAME\t$total_mapped_reads"
    
    for ref_list in $(ls ${REF_LIST_DIR}/*_contigs.txt); do
        ref_name=$(basename "$ref_list" _contigs.txt)
        # Robust column-1 match using awk (safer than grep -w when contig names contain dots/pipes)
        group_mapped_reads=$(awk 'NR==FNR{want[$1]=1; next} ($1 in want){s+=$3} END{print s+0}' \
            "$ref_list" "$IDXSTATS_FILE")
        if [ "$total_mapped_reads" -gt 0 ]; then
            proportion=$(awk -v g="$group_mapped_reads" -v t="$total_mapped_reads" 'BEGIN{printf "%.2f", (g/t)*100}')
        else
            proportion="0.00"
        fi
        proportions_line+="\t$ref_name:$group_mapped_reads($proportion%)"
    done
    
    echo -e "$proportions_line" >> "$PROPORTIONS_SUMMARY_FILE"
    echo "  - Proportions calculated and added to summary file for $BARCODE_NAME."

    # Clean up the temporary BED files
    rm -f "$TEMP_PLANT_BED" "$TEMP_SMAG_BED"

    # --- Task 3: Create clean SMAG FASTQ, including unmapped reads ---
    echo "  - Starting Task 3: Creating a clean FASTQ of SMAG and unmapped reads..."
    
    # Get the list of reads to keep
    READS_TO_KEEP_LIST="${SPLIT_OUTPUT_DIR}/smag_reads_to_keep.txt"
    
    # 1. Start with the list of reads that mapped uniquely to SMAG
    cat "${SPLIT_OUTPUT_DIR}/smag_unique_names.txt" > "$READS_TO_KEEP_LIST"

    # 2. Add the unmapped reads from the original BAM file
    # samtools view -f 4 filters for reads where the 'unmapped' flag is set
    samtools view -f 4 "$INPUT_BAM" | cut -f 1 >> "$READS_TO_KEEP_LIST"

    # 3. Remove duplicate entries
    sort -u "$READS_TO_KEEP_LIST" -o "$READS_TO_KEEP_LIST"

    # Filter the original FASTQ using the combined list
    INPUT_FASTQ="${FASTQ_DIR}/${BARCODE_NAME}.fastq.gz"
    CLEAN_FASTQ="${FINAL_OUTPUT_DIR}/${BARCODE_NAME}_smag_and_unmapped_clean.fastq.gz"

    if [ ! -f "$READS_TO_KEEP_LIST" ]; then
        echo "Error: Read name list not found for $BARCODE_NAME. Cannot proceed with Task 3."
        continue
    fi
    
    echo "  - Filtering reads from $INPUT_FASTQ..."
    seqtk subseq "$INPUT_FASTQ" "$READS_TO_KEEP_LIST" | gzip > "$CLEAN_FASTQ"
    echo "  - Task 3 finished for $BARCODE_NAME. Clean FASTQ saved to $CLEAN_FASTQ"

done

echo "All barcode processing complete. Summary is in $PROPORTIONS_SUMMARY_FILE"
echo "Job completed at $(date) on host $(hostname)"