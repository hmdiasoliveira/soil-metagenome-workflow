#!/bin/bash
#PBS -P <PROJECT>
#PBS -q gpuvolta
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/<PROJECT>+scratch/<PROJECT>
#PBS -l mem=384G
#PBS -l jobfs=400G
#PBS -l other=hyperthread
#PBS -l wd

module load cuda/12.8.0
set -euo pipefail
set -x

# ============================================================================
# Dorado Basecalling with Inline Demultiplexing
# ============================================================================
# Description: Basecall POD5 files using Dorado with adapter trimming and 
#              inline demultiplexing
# Input: POD5 files (directory)
# Output: Single BAM file with all barcodes
# ============================================================================

# Dorado paths
doradoBase="<DORADO_INSTALLATION_PATH>"
export LD_LIBRARY_PATH="${doradoBase}/lib:$LD_LIBRARY_PATH"
doradoBin="${doradoBase}/bin/dorado"

# Pinned model (directory form)
modelCfg="<DORADO_MODEL_PATH>"
export DORADO_MODELS_DIRECTORY="$(dirname "$modelCfg")"   # skip re-download

# I/O paths
InputDir="<INPUT_POD5_DIR>"
OutputDir="<OUTPUT_BAM_DIR>"
mkdir -p "$OutputDir"

OutputBam="${OutputDir}/$(basename "$InputDir")_trimmed.bam"

barcodeKit="<BARCODE_KIT>"

# Base-calling-
"$doradoBin" basecaller "$modelCfg" "$InputDir" \
    --device cuda:all \
    --kit-name "$barcodeKit" \
    --modified-bases 5mCG_5hmCG \
    --trim all \
    --recursive \
    > "$OutputBam"

echo "Job completed successfully at $(date)"