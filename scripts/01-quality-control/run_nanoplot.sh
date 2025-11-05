#!/bin/bash
#
# run_nanoplot.sh - Generate comprehensive read statistics with NanoPlot
#
# Usage:
#   ./run_nanoplot.sh <input_fastq> <output_dir> [threads]
#
# Arguments:
#   input_fastq : Path to input FASTQ file (can be .fastq or .fastq.gz)
#   output_dir  : Directory for NanoPlot output
#   threads     : Number of threads (default: 4)
#
# Example:
#   ./run_nanoplot.sh reads.fastq.gz qc/nanoplot 8
#

set -euo pipefail

# Check arguments
if [ "$#" -lt 2 ]; then
    echo "Error: Insufficient arguments"
    echo "Usage: $0 <input_fastq> <output_dir> [threads]"
    exit 1
fi

INPUT_FASTQ="$1"
OUTPUT_DIR="$2"
THREADS="${3:-4}"

# Validate input
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Error: Input file not found: $INPUT_FASTQ"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Extract sample name
SAMPLE_NAME=$(basename "$INPUT_FASTQ" .fastq.gz)
SAMPLE_NAME=${SAMPLE_NAME%.fastq}

echo "=================================================="
echo "Running NanoPlot"
echo "=================================================="
echo "Input:   $INPUT_FASTQ"
echo "Output:  $OUTPUT_DIR"
echo "Sample:  $SAMPLE_NAME"
echo "Threads: $THREADS"
echo "=================================================="

# Check if file has reads (skip empty files)
if [[ "$INPUT_FASTQ" == *.gz ]]; then
    HAS_READS=$(zcat "$INPUT_FASTQ" | head -n 1 | grep -c '^@' || true)
else
    HAS_READS=$(head -n 1 "$INPUT_FASTQ" | grep -c '^@' || true)
fi

if [ "$HAS_READS" -eq 0 ]; then
    echo "Warning: No reads found in $INPUT_FASTQ. Skipping."
    exit 0
fi

# Run NanoPlot
NanoPlot \
    --fastq "$INPUT_FASTQ" \
    --outdir "$OUTPUT_DIR" \
    --prefix "${SAMPLE_NAME}_" \
    --threads "$THREADS" \
    --loglength \
    --N50 \
    --plots hex dot kde \
    --format png pdf \
    --title "QC Report: $SAMPLE_NAME"

# Rename output files for clarity
if [ -f "${OUTPUT_DIR}/${SAMPLE_NAME}_NanoStats.txt" ]; then
    mv "${OUTPUT_DIR}/${SAMPLE_NAME}_NanoStats.txt" \
       "${OUTPUT_DIR}/${SAMPLE_NAME}.nanostats.txt"
fi

if [ -f "${OUTPUT_DIR}/${SAMPLE_NAME}_NanoPlot-report.html" ]; then
    mv "${OUTPUT_DIR}/${SAMPLE_NAME}_NanoPlot-report.html" \
       "${OUTPUT_DIR}/${SAMPLE_NAME}.nanoplot.html"
fi

echo ""
echo "NanoPlot completed successfully!"
echo "Results: $OUTPUT_DIR"