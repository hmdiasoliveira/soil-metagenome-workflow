#!/bin/bash

# --- Configuration ---
# Path to your Dorado v1.1.0 installation.
doradoBase="<DORADO_INSTALLATION_PATH>"
doradoBin="${doradoBase}/bin/dorado"

# The directory where you want to store the models.
models_dir="<DORADO_MODEL_DIR>"
mkdir -p "$models_dir"

# The specific model name to download.
# > simplex models
# - dna_r10.4.1_e8.2_400bps_hac@v5.2.0
# - dna_r10.4.1_e8.2_400bps_sup@v5.2.0
# ...
# > modification models
# - dna_r10.4.1_e8.2_400bps_hac@v5.2.0_6mA@v1
# - dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2
# ...

model_name="<DORADO_MODEL_NAME>"

# --- Execution ---
echo "Starting model download..."
echo "Downloading model: ${model_name} to directory: ${models_dir}"

# Use the dorado download command
"${doradoBin}" download --model "${model_name}" --models-directory "${models_dir}"

echo "Download script finished."