#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=5:00:00
#PBS -l mem=100GB
#PBS -l jobfs=100GB
#PBS -l ncpus=5
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Build Count Matrix from BAM Files
# ============================================================================
# Description: Extract read-to-contig mappings and build count matrix
# Input: Sorted BAM files from map_reads_to_assembly.sh
# Output: CSV with read mappings and count matrix
# ============================================================================

# Load samtools
module load samtools

# Input/Output paths
BAM_DIR="<INPUT_BAM_DIR>"
OUTPUT_CSV="<OUTPUT_MAPPED_READS_CSV>"
MATRIX_OUTPUT="<OUTPUT_COUNT_MATRIX_CSV>"

# Mapping quality threshold
MIN_MAPQ=60

# Create output directory
mkdir -p "$(dirname "$OUTPUT_CSV")"

echo "Starting count matrix generation..."

# Remove existing output file
if [ -f "$OUTPUT_CSV" ]; then
    rm "$OUTPUT_CSV"
fi

# Write CSV header
printf "sample,read_ID,contig_id,MAPQ\n" > "$OUTPUT_CSV"

echo "Processing BAM files..."

# Process each BAM file
for bam_file in "${BAM_DIR}"/*.sorted.bam; do
    if [ ! -f "$bam_file" ]; then
        echo "No BAM files found in ${BAM_DIR}"
        exit 1
    fi
    
    sample_name=$(basename "$bam_file" .sorted.bam)
    
    echo "Processing: ${sample_name}"
    
    # Extract mappings from BAM
    # -F 0x904: Filter unmapped, secondary, and supplementary alignments
    # -q: Minimum MAPQ threshold
    samtools view -@ "$PBS_NCPUS" -h -q "$MIN_MAPQ" -F 0x904 "$bam_file" \
        | awk -v sample_id="$sample_name" '
        BEGIN { FS="\t"; OFS="," }
        /^@/ { next }  # Skip header lines
        {
            read_id = $1
            contig_id = $3
            mapq = $5
            printf "%s,%s,%s,%s\n", sample_id, read_id, contig_id, mapq
        }' >> "$OUTPUT_CSV"
done

echo "Mapped reads table complete: $OUTPUT_CSV"

# Generate count matrix using R
echo "Generating count matrix..."

Rscript - <<EOF
# Read mapping data
in_file <- "$OUTPUT_CSV"
out_file <- "$MATRIX_OUTPUT"

cat("Reading mapped reads table...\n")

# Check if input exists
if (!file.exists(in_file)) {
    stop(paste("Error: Input file not found at", in_file))
}

# Read CSV
df <- read.csv(in_file, header = TRUE, stringsAsFactors = FALSE)

# Validate columns
if (!"sample" %in% names(df) || !"contig_id" %in% names(df)) {
    stop("Error: Required columns not found in CSV")
}

# Generate count matrix (samples × contigs)
cat("Building count matrix...\n")
count_matrix <- xtabs(~ sample + contig_id, data = df)

# Save matrix
write.csv(count_matrix, out_file, quote = FALSE)

cat("Count matrix saved to:", out_file, "\n")
cat(paste("Dimensions:", nrow(count_matrix), "samples ×", 
          ncol(count_matrix), "contigs\n"))
EOF

echo "Count matrix generation complete!"