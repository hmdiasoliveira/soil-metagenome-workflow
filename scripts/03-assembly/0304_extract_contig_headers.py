#!/usr/bin/env python3
"""
Extract Contig Headers from metaMDBG Assembly

Description: Parse metaMDBG contig headers and extract metadata
Input: Compressed FASTA file (contigs.fasta.gz)
Output: CSV file with contig metadata

Usage: python extract_contig_headers.py <input_fasta.gz> <output.csv>
"""

import gzip
import csv
import re
import sys

def extract_contig_headers(input_path, output_path):
    """
    Extract contig metadata from metaMDBG FASTA headers.
    
    Header format: >contig_id length=12345 coverage=67.89 circular=yes/no
    """
    # Regular expression to match header lines
    header_pattern = re.compile(r'^>(\S+)\s+length=(\d+)\s+coverage=([\d.]+)\s+circular=(\w+)')
    
    with gzip.open(input_path, 'rt') as infile, open(output_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["contig_id", "length", "coverage", "circular"])
        
        for line in infile:
            if line.startswith('>'):
                match = header_pattern.match(line.strip())
                if match:
                    writer.writerow(match.groups())

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_contig_headers.py <input_fasta.gz> <output.csv>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    print(f"Extracting headers from: {input_path}")
    extract_contig_headers(input_path, output_path)
    print(f"Headers saved to: {output_path}")