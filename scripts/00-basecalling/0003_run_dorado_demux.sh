#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=64GB
#PBS -l jobfs=100GB
#PBS -l ncpus=8
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Dorado Demultiplexing to FASTQ
# ============================================================================
# Description: Demultiplex basecalled BAM to per-barcode FASTQ files
# Input: Basecalled BAM file
# Output: Per-barcode FASTQ files (gzipped)
#
# IMPORTANT: Choose ONE of the two approaches below:
# 1. --no-classify: Use if BAM has existing barcode tags (BC:Z:)
# 2. --kit-name: Use to re-classify barcodes from scratch
# ============================================================================

# Load required modules
module load samtools/1.22
module load parallel
module load cuda/12.8.0
nvidia-smi || true

doradoBase="<DORADO_INSTALLATION_PATH>"
export LD_LIBRARY_PATH="${doradoBase}/lib:$LD_LIBRARY_PATH"
doradoBin="${doradoBase}/bin/dorado"

inputBAM="<INPUT_BAM_PATH>"
fastqDir="<OUTPUT_FASTQ_OR_BAM_DIR>"
barcodeKit="<BARCODE_KIT>"  # Only needed for Approach 2

mkdir -p "$fastqDir"

echo "Starting dorado demux..."

# ===== APPROACH 1: Use existing barcode tags (RECOMMENDED) =====
# (already basecalled and trimmed)
# Uncomment this if your BAM already has BC:Z: tags from basecalling
# Delete --emit-fastq if you only want BAM output with specific tags

"$doradoBin" demux \
  --output-dir "$fastqDir" \
  --no-classify \
  --emit-fastq \
  "$inputBAM"

# ===== APPROACH 2: Re-classify barcodes from scratch =====
# Uncomment this if barcode tags are missing or incorrect
# "$doradoBin" demux \
#   --kit-name "$barcodeKit" \
#   --output-dir "$fastqDir" \
#   --emit-fastq \
#   "$inputBAM"

echo "Compressing FASTQ files..."
find "$fastqDir" -name "*.fastq" -type f | parallel -j 12 'gzip {}'

echo "Complete! Gzipped FASTQ files are in $fastqDir"