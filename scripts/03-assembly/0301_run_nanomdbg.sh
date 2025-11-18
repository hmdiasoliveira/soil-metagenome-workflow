#!/bin/bash
#PBS -P <PROJECT>
#PBS -q hugemem
#PBS -l walltime=48:00:00
#PBS -l mem=1470GB
#PBS -l jobfs=1400GB
#PBS -l ncpus=48
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# metaMDBG Assembly - Nanopore Metagenome Assembly
# ============================================================================
# Description: Assemble Nanopore reads using metaMDBG (nanoMDBG method)
# Input: Filtered FASTQ file (quality and length filtered)
# Output: Assembled contigs in FASTA format
# ============================================================================

echo "Job started on: $(date)"
echo "Host: $(hostname)"
echo "PBS_JOBID: $PBS_JOBID"
echo "PBS_JOBFS: $PBS_JOBFS"

# Activate conda environment
echo "Activating conda environment..."
source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate <CONDA_ASSEMBLY_ENV> || { echo "ERROR: Conda activation failed!" >&2; exit 1; }

# Add metaMDBG and minimap2 to PATH
export PATH=<METAMDBG_PATH>/build/bin:$PATH
export PATH=<MINIMAP2_PATH>:$PATH

# Input and output paths
INPUT_READS="<INPUT_FILTERED_FASTQ>"
OUTPUT_DIR="<OUTPUT_ASSEMBLY_DIR>"
OUTPUT_BASENAME="metaMDBG_assembly"

SAMPLE_NAME=$(basename "$INPUT_READS" .fastq.gz)

echo "Processing: $SAMPLE_NAME"
echo "Input Reads: $INPUT_READS"
echo "Output Directory: $OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"

# Copy input to jobfs for faster I/O
LOCAL_INPUT="$PBS_JOBFS/$(basename "$INPUT_READS")"
LOCAL_OUTPUT_DIR="$PBS_JOBFS/metaMDBG_output"
mkdir -p "$LOCAL_OUTPUT_DIR"

echo "Copying input to jobfs..."
cp "$INPUT_READS" "$LOCAL_INPUT"

echo "Running metaMDBG assembly..."

# metaMDBG assembly with nanoMDBG method for Nanopore R10.4+ reads
#
# Basic options:
#   --out-dir               Output dir for contigs and temporary files
#   --in-ont                Nanopore R10.4+ read filename(s) (separated by space)
#   --threads               Number of cores
#
# Assembly options:
#   --kmer-size             k-mer size [15]
#   --density-assembly      Fraction of total k-mers used for assembly [0.005]
#   --max-k                 Stop assembly after k iterations [0]
#   --min-abundance         Minimum abundance for k-min-mers [0]
#                           Filter out unique k-min-mers to improve performance
#
# Correction options:
#   --min-read-quality      Minimum read average quality [0]
#   --density-correction    Fraction of total k-mers used for correction [0.025]
#   --min-read-identity     Min read identity [0.96]
#   --min-read-overlap      Min read overlap length [1000]

metaMDBG asm \
    --out-dir "$LOCAL_OUTPUT_DIR" \
    --in-ont "$LOCAL_INPUT" \
    --min-read-overlap 300 \
    --min-abundance 1 \
    --threads "$PBS_NCPUS" || { echo "ERROR: metaMDBG failed!" >&2; exit 1; }

# Capture environment and tool versions
echo "Capturing environment information..."
conda list > "$LOCAL_OUTPUT_DIR/conda_env_packages.txt"
metaMDBG --version > "$LOCAL_OUTPUT_DIR/metaMDBG_version.txt"

# Copy results to permanent storage
echo "Copying output to permanent storage..."
cp -r "$LOCAL_OUTPUT_DIR/"* "$OUTPUT_DIR/" || { echo "ERROR: Failed to copy results!" >&2; exit 1; }

echo "Assembly completed successfully at $(date)"