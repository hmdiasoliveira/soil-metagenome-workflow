#!/bin/bash
#PBS -P <PROJECT>
#PBS -q normal
#PBS -l walltime=05:00:00
#PBS -l mem=32GB
#PBS -l ncpus=8
#PBS -l storage=gdata/<PROJECT>+scratch/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Split Mixed FASTQ by Barcode Tags
# ============================================================================
# Description: Parse mixed FASTQ and separate reads by barcode from RG tags
#              Use this when dorado demux creates a single mixed FASTQ file
#              instead of per-barcode files
# Input: FASTQ file(s) with RG tags containing barcode information
# Output: Per-barcode FASTQ files (gzipped)
# ============================================================================

# Input mixed FASTQ file (from dorado demux output)
inputFASTQ="<INPUT_MIXED_FASTQ>"
outputDir="<OUTPUT_FASTQ_DIR>"

mkdir -p "$outputDir"

# Python script to split FASTQ by barcode
cat > /tmp/split_fastq_by_barcode.py << 'EOF'
#!/usr/bin/env python3
"""
Split a mixed FASTQ file into per-barcode files based on RG tags.

Reads FASTQ headers, extracts barcode information from the RG tag,
and writes each read to the appropriate barcode-specific output file.
"""

import sys
import gzip
import re
from collections import defaultdict

input_file = sys.argv[1]
output_dir = sys.argv[2]

# Dictionary to hold file handles for each barcode
barcode_files = {}
read_counts = defaultdict(int)

def get_barcode_from_header(header):
    """Extract barcode from RG tag in FASTQ header."""
    match = re.search(r'barcode(\d+)', header)
    if match:
        # Zero-pad barcode number to 2 digits (barcode01, barcode02, etc.)
        return f"barcode{match.group(1).zfill(2)}"
    else:
        return "unclassified"

# Open input file (handle both compressed and uncompressed)
if input_file.endswith('.gz'):
    infile = gzip.open(input_file, 'rt')
else:
    infile = open(input_file, 'r')

print(f"Processing: {input_file}", file=sys.stderr)

# Read FASTQ in chunks of 4 lines
line_count = 0
current_record = []

for line in infile:
    current_record.append(line)
    line_count += 1
    
    if line_count == 4:
        # Complete FASTQ record
        header = current_record[0]
        barcode = get_barcode_from_header(header)
        read_counts[barcode] += 1
        
        # Open output file for this barcode if not already open
        if barcode not in barcode_files:
            output_file = f"{output_dir}/{barcode}.fastq.gz"
            barcode_files[barcode] = gzip.open(output_file, 'wt')
            print(f"Created output file: {output_file}", file=sys.stderr)
        
        # Write record to appropriate file
        barcode_files[barcode].write(''.join(current_record))
        
        # Reset for next record
        current_record = []
        line_count = 0

infile.close()

# Close all output files
for f in barcode_files.values():
    f.close()

print(f"\nSplit into {len(barcode_files)} barcode files", file=sys.stderr)
print("\nRead counts by barcode:", file=sys.stderr)
for barcode in sorted(read_counts.keys()):
    print(f"  {barcode}: {read_counts[barcode]:,} reads", file=sys.stderr)

EOF

# Run the Python script
echo "Splitting FASTQ by barcode..."
python3 /tmp/split_fastq_by_barcode.py "$inputFASTQ" "$outputDir"

# Report results
echo ""
echo "============================================"
echo "Demultiplexing complete!"
echo "============================================"
echo ""
echo "Output directory: $outputDir"
echo ""
echo "Files created:"
ls -lh "$outputDir"/*.fastq.gz

echo ""
echo "Read counts per barcode:"
for fq in "$outputDir"/*.fastq.gz; do
    count=$(zcat "$fq" | awk 'NR%4==1' | wc -l)
    printf "  %-30s %10s reads\n" "$(basename $fq):" "$count"
done

echo ""
total=$(zcat "$outputDir"/*.fastq.gz | awk 'NR%4==1' | wc -l)
echo "Total reads: $total"

# Clean up temporary Python script
rm -f /tmp/split_fastq_by_barcode.py