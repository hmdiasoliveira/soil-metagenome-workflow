#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=64GB
#PBS -l ncpus=16
#PBS -l storage=gdata/<PROJECT>+scratch/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Split Barcoded BAM by Read Groups (Alternative to dorado demux)
# ============================================================================
# Description: Separate inline-barcoded BAM into per-barcode FASTQ files
#              using samtools split. Use this when dorado demux creates
#              nested directories instead of per-barcode files.
# Input: BAM file with RG tags containing barcode information
# Output: Per-barcode FASTQ files (gzipped)
# ============================================================================

module load samtools/1.19
module load parallel

# Input BAM with inline barcode tags
inputBAM="<INPUT_BARCODED_BAM>"
tempBAMDir="<TEMP_BAM_DIR>"
fastqDir="<OUTPUT_FASTQ_DIR>"

mkdir -p "$tempBAMDir"
mkdir -p "$fastqDir"

# Step 1: Split BAM by read group (which contains barcode info)
echo "Splitting BAM by read groups..."
samtools split \
    -@ 16 \
    -u "$tempBAMDir/unclassified.bam" \
    -f "$tempBAMDir/RG_%!.bam" \
    "$inputBAM"

# Step 2: Convert each BAM to FASTQ and rename properly
echo "Converting BAMs to FASTQ and extracting barcode numbers..."

for bamFile in "$tempBAMDir"/RG_*.bam; do
    # Get the actual RG string from the BAM header
    rgString=$(basename "$bamFile" .bam | sed 's/RG_//')
    
    # Extract barcode number from RG string using the BAM header
    barcode=$(samtools view -H "$bamFile" | grep '@RG' | grep -oP 'barcode\d+' | head -n1)
    
    if [ -z "$barcode" ]; then
        # Fallback: try to extract from RG string itself
        barcode=$(echo "$rgString" | grep -oP 'barcode\d+' || echo "unknown_$(basename $bamFile .bam)")
    fi
    
    echo "  Processing: $rgString -> $barcode"
    
    # Convert BAM to FASTQ and compress
    samtools fastq -@ 2 -n "$bamFile" | gzip -c > "$fastqDir/${barcode}.fastq.gz"
done

# Step 3: Handle unclassified reads
if [ -f "$tempBAMDir/unclassified.bam" ]; then
    echo "  Processing unclassified reads..."
    samtools fastq -@ 2 -n "$tempBAMDir/unclassified.bam" | gzip -c > "$fastqDir/unclassified.fastq.gz"
fi

# Step 4: Merge duplicate barcode files if multiple RGs have same barcode
echo "Checking for duplicate barcodes from different read groups..."
cd "$fastqDir"

for barcode in $(ls *.fastq.gz | sed 's/\.fastq\.gz$//' | sort -u); do
    # Count files matching this barcode pattern
    matching_files=$(ls ${barcode}*.fastq.gz 2>/dev/null | wc -l)
    
    if [ "$matching_files" -gt 1 ]; then
        echo "  Merging $matching_files files for $barcode"
        zcat ${barcode}*.fastq.gz | gzip > ${barcode}_merged.fastq.gz
        rm ${barcode}*.fastq.gz
        mv ${barcode}_merged.fastq.gz ${barcode}.fastq.gz
    fi
done

# Step 5: Clean up temporary BAM files
echo "Cleaning up temporary BAM files..."
rm -rf "$tempBAMDir"

# Step 6: Report results
echo ""
echo "============================================"
echo "Demultiplexing complete!"
echo "============================================"
echo ""
echo "Output directory: $fastqDir"
echo ""
echo "Files created:"
ls -lh "$fastqDir"/*.fastq.gz

echo ""
echo "Read counts per barcode:"
for fq in "$fastqDir"/*.fastq.gz; do
    count=$(zcat "$fq" | awk 'NR%4==1' | wc -l)
    printf "  %-30s %10s reads\n" "$(basename $fq):" "$count"
done

echo ""
echo "Total reads:"
total=$(zcat "$fastqDir"/*.fastq.gz | awk 'NR%4==1' | wc -l)
echo "  $total"