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
# metaMDBG Polishing - Polish Assembled Contigs
# ============================================================================
# Description: Polish metaMDBG contigs using original Nanopore reads
# Input: Assembled contigs (FASTA) and original reads (FASTQ)
# Output: Polished contigs in FASTA format
#
# NOTE: The standalone 'metaMDBG polish' subcommand is only available in
#       metaMDBG v1.2. In v1.3+, polishing is integrated into the 'asm'
#       command and contigs are automatically polished during assembly.
#       If you are using v1.3+, this script is not needed — use
#       run_nanomdbg.sh instead.
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
export PATH=<METAMDBG_V1.2_PATH>/build/bin:$PATH
export PATH=<MINIMAP2_PATH>:$PATH

# Input files
INPUT_CONTIGS="<INPUT_CONTIGS_FASTA>"
INPUT_READS="<INPUT_FILTERED_FASTQ>"

# Output directory
OUTPUT_DIR="<OUTPUT_POLISH_DIR>"
mkdir -p "$OUTPUT_DIR"

echo "Input contigs: $INPUT_CONTIGS"
echo "Input reads: $INPUT_READS"
echo "Output directory: $OUTPUT_DIR"

# metaMDBG polish (v1.2 only)
#
# Basic options:
#   --out-dir               Output dir for polished contigs
#   --polish-target         Contigs to polish
#   --in-ont                Nanopore R10.4+ read filename(s)
#   --threads               Number of cores
#
# Other options:
#   -n                      Maximum read coverage for correction [0]
#   --max-memory            Maximum memory usage for read mapping [8]

metaMDBG polish \
    --out-dir "$OUTPUT_DIR" \
    --polish-target "$INPUT_CONTIGS" \
    --in-ont "$INPUT_READS" \
    --threads "$PBS_NCPUS" || { echo "ERROR: Polishing failed!" >&2; exit 1; }

echo "Polishing completed successfully at $(date)"