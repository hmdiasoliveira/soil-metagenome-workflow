#!/bin/bash
#PBS -P <PROJECT>
#PBS -q dgxa100
#PBS -l ncpus=16
#PBS -l ngpus=1
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/<PROJECT>+scratch/<PROJECT>
#PBS -l mem=384G
#PBS -l jobfs=1200G
#PBS -l other=hyperthread
#PBS -l wd

module load cuda/12.8.0
set -euo pipefail
set -x

# ============================================================================
# Resume Dorado Basecalling from Incomplete Run
# ============================================================================
# Description: Resume basecalling from an incomplete BAM file
# Input: POD5 files (directory) + incomplete BAM file
# Output: Completed BAM file with additional basecalled reads
# ============================================================================

# --- Dorado paths -------------------------------------------------
doradoBase="<DORADO_INSTALLATION_PATH>"
export LD_LIBRARY_PATH="${doradoBase}/lib:$LD_LIBRARY_PATH"
doradoBin="${doradoBase}/bin/dorado"

# Pinned model (directory form)
modelCfg="<DORADO_MODEL_PATH>"
export DORADO_MODELS_DIRECTORY="$(dirname "$modelCfg")"

# --- I/O paths and file names -------------------------------------
InputDir="<INPUT_POD5_DIR>"
OutputDir="<OUTPUT_BAM_DIR>"

# Define the paths to the incomplete and new resumed BAM files
incompleteBam="<INCOMPLETE_BAM_PATH>"
resumedBam="<RESUMED_BAM_OUTPUT_PATH>"

barcodeKit="<BARCODE_KIT>"

# --- Validation ---------------------------------------------------
# Check if the file to resume from exists
if [ ! -f "$incompleteBam" ]; then
    echo "Error: Incomplete BAM file not found at $incompleteBam"
    exit 1
fi

echo "Resuming basecalling from incomplete BAM file: $incompleteBam"

# --- Base-calling -------------------------------------------------
# Resume base-calling using --resume-from and output to the new file
"$doradoBin" basecaller "$modelCfg" "$InputDir" \
    --device cuda:all \
    --kit-name "$barcodeKit" \
    --trim all \
    --recursive \
    --resume-from "$incompleteBam" \
    > "$resumedBam"

echo "Job completed successfully at $(date)"