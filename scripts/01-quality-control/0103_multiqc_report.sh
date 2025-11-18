#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# MultiQC Report Generation
# ============================================================================
# Description: Aggregate FastQC and NanoPlot results into single report
# Input: FastQC and NanoPlot output directories
# Output: MultiQC HTML report
# ============================================================================

# Activate conda environment
source <CONDA_PATH>/etc/profile.d/conda.sh
conda activate <QC_ENV_PATH>

# Input/Output paths (MODIFY THESE)
QC_BASE="<QC_RESULTS_BASE_DIR>"
FASTQC_DIR="${QC_BASE}/fastqc"
NANOPLOT_DIR="${QC_BASE}/nanoplot"
OUT_MULTIQC="${QC_BASE}/multiqc"

mkdir -p "$OUT_MULTIQC"

# Check if input directories exist
if [ ! -d "$FASTQC_DIR" ] && [ ! -d "$NANOPLOT_DIR" ]; then
    echo "Error: Neither FastQC nor NanoPlot directories found"
    echo "  FastQC: $FASTQC_DIR"
    echo "  NanoPlot: $NANOPLOT_DIR"
    exit 1
fi

echo "Generating MultiQC report..."
echo "  FastQC results: $FASTQC_DIR"
echo "  NanoPlot results: $NANOPLOT_DIR"
echo "  Output: $OUT_MULTIQC"

# Run MultiQC with original options
multiqc \
    "$FASTQC_DIR" \
    "$NANOPLOT_DIR" \
    -o "$OUT_MULTIQC"

echo "MultiQC report generated: ${OUT_MULTIQC}/multiqc_report.html"
echo "MultiQC completed at $(date) on host $(hostname)"