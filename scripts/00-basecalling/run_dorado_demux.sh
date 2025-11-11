#!/bin/bash
#PBS -P <PROJECT>
#PBS -q dgxa100
#PBS -l ncpus=16
#PBS -l ngpus=1
#PBS -l walltime=20:00:00
#PBS -l storage=gdata/<PROJECT>+scratch/<PROJECT>
#PBS -l mem=96G
#PBS -l jobfs=100G
#PBS -l other=hyperthread
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Dorado Demultiplexing and BAM to FASTQ Conversion
# ============================================================================
# Description: Demultiplex basecalled BAM and convert to gzipped FASTQ files
# Input: Basecalled BAM file (single file or directory)
# Output: Per-barcode BAM files and gzipped FASTQ files
# ============================================================================

# Load required modules
module load samtools/1.19
module load parallel
module load cuda/12.8.0
nvidia-smi || true

# Paths to dorado binaries and environment setup
doradoBase="<DORADO_INSTALLATION_PATH>"
export LD_LIBRARY_PATH="${doradoBase}/lib:$LD_LIBRARY_PATH"
doradoBin="${doradoBase}/bin/dorado"

# Input BAM file (already basecalled and trimmed)
inputBAM="<INPUT_BAM_PATH>"
demuxDir="<OUTPUT_DEMUX_BAM_DIR>"
fastqDir="<OUTPUT_FASTQ_DIR>"

# Create output directories
mkdir -p "$demuxDir"
mkdir -p "$fastqDir"

echo "Starting dorado demux..."
"$doradoBin" demux \
  --output-dir "$demuxDir" \
  --no-classify \
  --emit-fastq \
  "$inputBAM"

echo "Demux complete. Converting BAM to gzipped FASTQ..."

# Convert each barcode BAM in output to gzipped FASTQ
find "$demuxDir" -name "*.bam" | parallel -j 12 "samtools fastq -@ 1 -n {} | gzip -c > $fastqDir/{/.}.fastq.gz"

echo "All processes complete. FASTQ files are in $fastqDir"