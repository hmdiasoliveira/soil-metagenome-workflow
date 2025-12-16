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

# Model configuration (LOCAL ONLY; no downloading in PBS job)

# Directory where all Dorado models are stored locally (downloaded beforehand)
modelsDir="<DORADO_MODEL_DIR>"

# Simplex model directory (basecalling A/C/G/T)
simplexModelPath="<DORADO_SIMPLEX_MODEL_PATH>"

# Modification model directory (e.g., 5mCG_5hmCG compatible with the simplex model)
modModelPath="<DORADO_MOD_MODEL_PATH>"

# Make Dorado look here for models (prevents auto-download attempts)
export DORADO_MODELS_DIRECTORY="$modelsDir"

# I/O paths
InputDir="<INPUT_POD5_DIR>"
OutputDir="<OUTPUT_BAM_DIR>"
mkdir -p "$OutputDir"

OutputBam="${OutputDir}/$(basename "$InputDir")_trimmed.bam"

barcodeKit="<BARCODE_KIT>"

# Base-calling (simplex + explicit modified-base model path)

"$doradoBin" basecaller "$simplexModelPath" "$InputDir" \
  --models-directory "$modelsDir" \
  --modified-bases-models "$modModelPath" \
  --device cuda:all \
  --kit-name "$barcodeKit" \
  --trim all \
  --recursive \
  > "$OutputBam"

echo "Job completed successfully at $(date)"