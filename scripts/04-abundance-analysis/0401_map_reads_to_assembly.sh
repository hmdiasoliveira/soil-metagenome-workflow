#!/bin/bash
#PBS -P <PROJECT>
#PBS -q hugemem
#PBS -l walltime=24:00:00
#PBS -l mem=500GB
#PBS -l jobfs=250GB
#PBS -l ncpus=48
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Map Reads to Assembly - Minimap2 Alignment
# ============================================================================
# Description: Map filtered reads back to assembled contigs for abundance
# Input: Assembled contigs (FASTA) and filtered reads (FASTQ)
# Output: Sorted and indexed BAM files per sample
# ============================================================================

# Add minimap2 to PATH
export PATH=<MINIMAP2_PATH>:$PATH

# Load samtools
module load samtools

# Input files
CONTIGS_FILE="<INPUT_CONTIGS_FASTA>"
READS_DIR="<INPUT_FILTERED_FASTQ_DIR>"
OUTPUT_DIR="<OUTPUT_BAM_DIR>"
THREADS=$PBS_NCPUS

# Create output directory
echo "Creating output directory: ${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Map reads for all samples
echo "Starting read mapping for all samples..."

for read_file in "${READS_DIR}"/*.fastq.gz; do
    # Get sample name
    base_name=$(basename "${read_file}" .fastq.gz)
    output_bam="${OUTPUT_DIR}/${base_name}.sorted.bam"
    
    echo "Processing sample: ${base_name}"
    echo "Output BAM: ${output_bam}"
    
    # Map reads with minimap2 and sort with samtools
    # -ax map-ont: Preset for Oxford Nanopore reads
    # -t: Number of threads
    # Pipe to samtools sort for immediate sorting
    minimap2 -ax map-ont -t "${THREADS}" "${CONTIGS_FILE}" "${read_file}" \
        | samtools sort -o "${output_bam}" -@ "${THREADS}" -
    
    # Check if sorting was successful
    if [ $? -eq 0 ]; then
        echo "Successfully mapped and sorted reads for ${base_name}"
        
        # Index the BAM file
        echo "Indexing BAM file for ${base_name}..."
        samtools index -c "${output_bam}"
        
        if [ $? -eq 0 ]; then
            echo "Successfully indexed BAM file for ${base_name}"
        else
            echo "ERROR: Failed to index BAM for ${base_name}"
            exit 1
        fi
    else
        echo "ERROR: Failed to map/sort reads for ${base_name}"
        exit 1
    fi
done

echo "All read mapping tasks complete!"