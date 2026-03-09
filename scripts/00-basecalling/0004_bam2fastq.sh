#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=64GB
#PBS -l jobfs=100GB
#PBS -l ncpus=8
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

# Load any required modules (none needed for dorado bam2fastq unless on path)
module load samtools/1.22
module load parallel

set -euo pipefail
set -x

# Directory containing per-barcode BAM files (e.g. barcode01.bam, barcode02.bam, ...)
inputDir="<INPUT_BAM_PATH>"
outputDir="<OUTPUT_FASTQ_DIR>"

mkdir -p "$outputDir"

find "$inputDir" -name "*.bam" | parallel -j 12 "samtools fastq -T MM,ML -@ 1 -n {} | gzip -c > $outputDir/{/.}.fastq.gz"

echo "FASTQ generation complete → files saved in $outputDir"