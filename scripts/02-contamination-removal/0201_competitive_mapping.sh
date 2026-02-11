#!/bin/bash
#PBS -P <PROJECT>
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l mem=1000GB
#PBS -l jobfs=500GB
#PBS -l ncpus=48
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Competitive Mapping to Combined Plant + SMAG Reference
# ============================================================================
# Description: Map reads to combined reference containing both plant genomes
#              and microbial genomes (SMAG) to identify host contamination
# Input: Filtered FASTQ files (gzipped)
# Output: Sorted BAM files with primary and secondary alignments
# ============================================================================

# Use custom minimap2 version (manual install)
export PATH=<MINIMAP2_PATH>:$PATH

# Load samtools module (for BAM processing)
module load samtools/1.19

# Reference genome prebuilt index (.mmi) and dict file for @SQ headers
REFERENCE_MMI="<COMBINED_REFERENCE_MMI>"
DICT_FILE="<COMBINED_REFERENCE_DICT>"

# Input FASTQ directory (gzipped reads)
INPUT_DIR="<INPUT_FASTQ_DIR>"

# Output BAM directory
OUTPUT_DIR="<OUTPUT_BAM_DIR>"
mkdir -p "$OUTPUT_DIR"

# Copy FASTQ and reference index to local jobfs for faster I/O
echo "Copying input files to jobfs..."
cp -r "$INPUT_DIR" "$PBS_JOBFS"
cp "$REFERENCE_MMI" "$PBS_JOBFS"

LOCAL_INPUT_DIR="$PBS_JOBFS/$(basename "$INPUT_DIR")"
LOCAL_REFERENCE_MMI="$PBS_JOBFS/$(basename "$REFERENCE_MMI")"

# Mapping loop
# Retain all primary mappings. Also retain up to -N [=5] top secondary mappings 
# if their chaining scores are higher than -p [=0.8] of their corresponding 
# primary mappings.
find "$LOCAL_INPUT_DIR" -name "*.fastq.gz" | while read -r FASTQ_FILE
do
    SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq.gz)
    SAM_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.sam"
    BAM_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.sorted.bam"

    echo "Processing $SAMPLE_NAME..."

    # Mapping → SAM output
    # Map ONT reads; keep secondaries for host filtering; include MD/NM tags.
    # --secondary=yes [default]
    # Minimap2 rates an alignment by the score of the max-scoring sub-segment, 
    # excluding introns, and marks the best alignment as primary in SAM.
    minimap2 \
        -t $PBS_NCPUS \
        -ax map-ont \
        --MD \
        "$LOCAL_REFERENCE_MMI" \
        "$FASTQ_FILE" \
        > "$SAM_FILE"

    # Check if @SQ headers exist
    if ! grep -q "^@SQ" "$SAM_FILE"; then
        echo "Adding @SQ headers to $SAM_FILE"
        cat "$DICT_FILE" "$SAM_FILE" > "${OUTPUT_DIR}/${SAMPLE_NAME}_fixed.sam"
        mv "${OUTPUT_DIR}/${SAMPLE_NAME}_fixed.sam" "$SAM_FILE"
    fi

    # Convert to BAM + sort
    samtools view -@ $PBS_NCPUS -bS "$SAM_FILE" |
    samtools sort -@ $PBS_NCPUS -o "$BAM_FILE"

    # Index BAM
    samtools index -c "$BAM_FILE"
    
    # Mapping summary
    samtools flagstat "$BAM_FILE" > "${BAM_FILE%.bam}.flagstat.txt"
    samtools stats "$BAM_FILE" > "${BAM_FILE%.bam}.stats.txt"
    samtools idxstats "$BAM_FILE" > "${BAM_FILE%.bam}.idxstats.txt"

    # Remove SAM to save space
    rm "$SAM_FILE"

    echo "Finished $SAMPLE_NAME → BAM and index ready."

done

echo "All samples processed → BAMs saved in $OUTPUT_DIR"
echo "Job completed at $(date) on host $(hostname)"
