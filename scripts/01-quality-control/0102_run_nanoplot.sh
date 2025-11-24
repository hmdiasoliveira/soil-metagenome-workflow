#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=05:00:00
#PBS -l mem=64GB
#PBS -l jobfs=400GB
#PBS -l ncpus=8
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# NanoPlot Quality Assessment
# ============================================================================
# Description: Generate comprehensive read statistics and quality visualizations
# Input: FASTQ files (gzipped or uncompressed)
# Output: NanoStats.txt and NanoPlot-report.html per sample
# ============================================================================

# Activate conda environment
source <CONDA_PATH>/etc/profile.d/conda.sh
conda activate <QC_ENV_PATH>

# Configuration
NTHREADS=${PBS_NCPUS:-8}

# Input/Output paths (MODIFY THESE)
INPUT_DIR="<INPUT_FASTQ_DIR>"
OUT_BASE="<OUTPUT_BASE_DIR>"
OUT_NANOPLOT="${OUT_BASE}/nanoplot"

mkdir -p "$OUT_NANOPLOT"

# Local temporary directory on jobfs (faster I/O)
LOCAL_DIR="${PBS_JOBFS}/nanoplot_work"
LOCAL_INPUT="${LOCAL_DIR}/input"
LOCAL_NANOPLOT="${LOCAL_DIR}/nanoplot_results"

mkdir -p "$LOCAL_INPUT" "$LOCAL_NANOPLOT"

# Copy input files to local jobfs
echo "Copying input files to jobfs..."
cp "$INPUT_DIR"/*.fastq.gz "$LOCAL_INPUT/" 2>/dev/null || \
cp "$INPUT_DIR"/*.fq.gz "$LOCAL_INPUT/" 2>/dev/null || \
cp "$INPUT_DIR"/*.fastq "$LOCAL_INPUT/" 2>/dev/null || \
cp "$INPUT_DIR"/*.fq "$LOCAL_INPUT/" 2>/dev/null || \
{ echo "No FASTQ files found in $INPUT_DIR"; exit 1; }

# Process each FASTQ file
for fq_file in "$LOCAL_INPUT"/*.f*q*; do
    # Extract sample name (remove extensions)
    sample=$(basename "$fq_file" .fastq.gz)
    sample=${sample%.fq.gz}
    sample=${sample%.fastq}
    sample=${sample%.fq}
    
    echo "Processing $sample..."
    
    # Check if file contains reads
    if gzip -dc "$fq_file" 2>/dev/null | grep -q '^@'; then
        
        # Create sample-specific output directory
        nanoplot_sample_dir="${LOCAL_NANOPLOT}/nanoplot_${sample}"
        mkdir -p "$nanoplot_sample_dir"
        
        # Run NanoPlot with original options
        NanoPlot \
            --fastq "$fq_file" \
            -o "$nanoplot_sample_dir" \
            --threads "$NTHREADS" \
            --loglength \
            --N50 \
            --plots dot
        
        # Move results to final output with sample-specific names
        if [ -f "${nanoplot_sample_dir}/NanoStats.txt" ]; then
            mv "${nanoplot_sample_dir}/NanoStats.txt" "${OUT_NANOPLOT}/${sample}.txt"
        fi
        
        if [ -f "${nanoplot_sample_dir}/NanoPlot-report.html" ]; then
            mv "${nanoplot_sample_dir}/NanoPlot-report.html" "${OUT_NANOPLOT}/${sample}.html"
        fi
        
        echo "Completed processing $sample"
    else
        echo "Skipping $sample: no reads found"
    fi
done

# Cleanup
rm -rf "$LOCAL_DIR"

echo "NanoPlot completed at $(date) on host $(hostname)"