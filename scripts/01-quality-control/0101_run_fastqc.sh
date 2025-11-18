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
# FastQC Quality Assessment
# ============================================================================
# Description: Per-base sequence quality analysis
# Input: FASTQ files (gzipped or uncompressed)
# Output: FastQC HTML reports and ZIP files per sample
# ============================================================================

# Activate conda environment
source <CONDA_PATH>/etc/profile.d/conda.sh
conda activate <QC_ENV_PATH>

# Configuration
NTHREADS=${PBS_NCPUS:-8}

# Increase Java heap space for FastQC
export _JAVA_OPTIONS="-Xmx48G"

# Input/Output paths (MODIFY THESE)
INPUT_DIR="<INPUT_FASTQ_DIR>"
OUT_BASE="<OUTPUT_BASE_DIR>"
OUT_FASTQC="${OUT_BASE}/fastqc"

mkdir -p "$OUT_FASTQC"

# Local temporary directory on jobfs (faster I/O)
LOCAL_DIR="${PBS_JOBFS}/fastqc_work"
LOCAL_INPUT="${LOCAL_DIR}/input"
LOCAL_FASTQC="${LOCAL_DIR}/fastqc_results"

mkdir -p "$LOCAL_INPUT" "$LOCAL_FASTQC"

# Copy input files to local jobfs
echo "Copying input files to jobfs..."
cp "$INPUT_DIR"/*.fastq.gz "$LOCAL_INPUT/" 2>/dev/null || \
cp "$INPUT_DIR"/*.fq.gz "$LOCAL_INPUT/" 2>/dev/null || \
cp "$INPUT_DIR"/*.fastq "$LOCAL_INPUT/" 2>/dev/null || \
cp "$INPUT_DIR"/*.fq "$LOCAL_INPUT/" 2>/dev/null || \
{ echo "No FASTQ files found in $INPUT_DIR"; exit 1; }

# Run FastQC on all files with original options
echo "Running FastQC..."
fastqc \
    "$LOCAL_INPUT"/*.f*q* \
    -o "$LOCAL_FASTQC" \
    -t "$NTHREADS"

# Copy results to final output
echo "Copying results to output directory..."
cp "$LOCAL_FASTQC"/* "$OUT_FASTQC"/

# Cleanup
rm -rf "$LOCAL_DIR"

echo "FastQC completed at $(date) on host $(hostname)"