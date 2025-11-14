#!/usr/bin/env python3
"""
Classify reads as plant-unique, SMAG-unique, or ambiguous based on competitive alignment.

This script analyzes BAM files from competitive mapping to a combined plant+SMAG reference.
For each read, it compares the best alignment score to plant contigs vs. SMAG contigs
and assigns the read to a category based on tie-breaking logic.

Usage:
    python classify_reads.py <input_bam> <ref_list_dir> <output_dir>

Arguments:
    input_bam:    BAM file from competitive mapping (sorted, indexed)
    ref_list_dir: Directory containing plant_contigs.txt and smag_contigs.txt
    output_dir:   Output directory for read name lists

Output:
    plant_unique_names.txt: Reads that map better to plant contigs
    smag_unique_names.txt:  Reads that map better to SMAG contigs
    ambiguous_names.txt:    Reads with equal scores to both categories
"""

import sys
import os
import pysam
from collections import defaultdict

# ----------------------------- Optional knobs -----------------------------
# Minimum MAPQ to consider an alignment (0 keeps everything).
MIN_MAPQ = 0

# Whether to include secondary / supplementary alignments.
# Your original script implicitly INCLUDED them (you didn't skip),
# which is reasonable for ONT multi-mapping. Keep that default here.
INCLUDE_SECONDARY = True
INCLUDE_SUPPLEMENTARY = True
# --------------------------------------------------------------------------


def _score_tuple(read):
    """
    Build a robust tie-break score for an alignment record.
    Order of comparison: AS (alignment score), MAPQ, aligned length.

    Returns a tuple (AS, MAPQ, aligned_length) where:
      - AS defaults to a very low value if absent to avoid artificially
        winning ties against real-scored alignments.
      - MAPQ defaults to 0 if missing.
      - aligned_length defaults to 0 if missing.
    """
    # Alignment score (AS) may be absent in some pipelines
    ascore = read.get_tag('AS') if read.has_tag('AS') else -10**9
    # Mapping quality (int, can be 0); None -> 0
    mapq = read.mapping_quality if read.mapping_quality is not None else 0
    # Query-aligned length (int; can be 0); None -> 0
    alen = read.query_alignment_length or 0
    return (ascore, mapq, alen)


def main(input_bam_path, ref_list_dir, output_dir):
    """
    Splits reads from a BAM file into plant-unique, SMAG-unique, and ambiguous
    based on a tie-breaking logic.
    
    The script iterates through all alignments in the BAM file. For each read,
    it keeps track of the best alignment score for any hit to a plant contig
    and the best score for any hit to a SMAG contig. After evaluating all alignments
    for all reads, it compares the best plant score against the best SMAG score
    for each read to make a final assignment.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # -------------------------------------------------------------------------
    # PART 1: Load reference contig names from text files
    # -------------------------------------------------------------------------
    # We load the lists of contig names for plants and SMAGs into sets.
    # Using sets for this lookup is very fast (O(1) average time complexity),
    # which is crucial for iterating over a large number of reads.
    plant_contigs = set()
    with open(f"{ref_list_dir}/plant_contigs.txt", 'r') as f:
        for line in f:
            name = line.strip()
            if name:
                plant_contigs.add(name)

    smag_contigs = set()
    with open(f"{ref_list_dir}/smag_contigs.txt", 'r') as f:
        for line in f:
            name = line.strip()
            if name:
                smag_contigs.add(name)
            
    # -------------------------------------------------------------------------
    # PART 2: Iterate through the BAM and find the best alignment for each read
    # -------------------------------------------------------------------------
    
    # Dictionaries to hold the best alignment quality for a given read
    # for each reference category (plant or SMAG). The key is the read name,
    # and the value is a tuple representing the best score found so far.
    # The score is a tuple (alignment_score, MAPQ, CIGAR_length) for a more robust
    # comparison than just a single value. Python compares tuples element by element.
    best_smag_alignments = {}
    best_plant_alignments = {}

    # This dictionary will store the final category for each unique read name.
    # We will populate this after we have processed all alignments for all reads.
    read_categories = {}

    # Counters for logging
    seen_alignments = 0
    skipped_unmapped = 0
    skipped_mapq = 0
    skipped_flags = 0

    # Open the input BAM file for reading. The 'rb' flag means "read binary".
    with pysam.AlignmentFile(input_bam_path, "rb") as bamfile:
        # Iterate over every alignment record in the BAM file.
        # `fetch(until_eof=True)` ensures we read the entire file.
        for read in bamfile.fetch(until_eof=True):
            seen_alignments += 1

            # Skip unmapped reads as they are not relevant for this analysis.
            if read.is_unmapped:
                skipped_unmapped += 1
                continue

            # Optional: filter by MAPQ
            if MIN_MAPQ and (read.mapping_quality is None or read.mapping_quality < MIN_MAPQ):
                skipped_mapq += 1
                continue

            # Optional: include/exclude secondary/supplementary alignments
            if (not INCLUDE_SECONDARY and read.is_secondary) or \
               (not INCLUDE_SUPPLEMENTARY and read.is_supplementary):
                skipped_flags += 1
                continue

            read_id = read.query_name

            # Build the comparison tuple
            score = _score_tuple(read)

            # Check which reference category this alignment belongs to.
            # This is where we use the sets created in Part 1 for quick lookups.
            ref_name = read.reference_name
            is_plant_hit = ref_name in plant_contigs
            is_smag_hit  = ref_name in smag_contigs
            
            # Update the best alignment score for this read in the plant category.
            # The logic checks if this is the first hit or if the current score
            # (alignment_score, MAPQ, cigar_length) is better than the previously recorded best.
            if is_plant_hit:
                prev = best_plant_alignments.get(read_id)
                if prev is None or score > prev:
                    best_plant_alignments[read_id] = score
            
            # Do the same for the SMAG category.
            if is_smag_hit:
                prev = best_smag_alignments.get(read_id)
                if prev is None or score > prev:
                    best_smag_alignments[read_id] = score
    
    # -------------------------------------------------------------------------
    # PART 3: Make the final decision and write read names to files
    # -------------------------------------------------------------------------

    # Now, iterate through all unique read IDs that had at least one hit.
    # The union of the keys from both dictionaries gives us all such read IDs.
    all_ids = set(best_plant_alignments.keys()) | set(best_smag_alignments.keys())

    # Open the three output files for writing the read names.
    plant_path = f"{output_dir}/plant_unique_names.txt"
    smag_path = f"{output_dir}/smag_unique_names.txt"
    ambi_path = f"{output_dir}/ambiguous_names.txt"

    plant_count = smag_count = ambi_count = 0

    with open(plant_path, 'w') as f_plant, \
         open(smag_path, 'w') as f_smag, \
         open(ambi_path, 'w') as f_ambiguous:

        for read_id in all_ids:
            # Retrieve the best score for each category. Use a very low tuple if absent
            # so that a real score always wins.
            plant_score = best_plant_alignments.get(read_id, (-10**9, 0, 0))
            smag_score  = best_smag_alignments.get(read_id,  (-10**9, 0, 0))
            
            # Apply the tie-breaking logic.
            # 1. If the read only had a plant hit (and it was a good one), it's plant-unique.
            if plant_score > (-10**9, 0, 0) and smag_score == (-10**9, 0, 0):
                f_plant.write(f"{read_id}\n")
                plant_count += 1
            # 2. If the read only had a SMAG hit, it's SMAG-unique.
            elif smag_score > (-10**9, 0, 0) and plant_score == (-10**9, 0, 0):
                f_smag.write(f"{read_id}\n")
                smag_count += 1
            # 3. If it had hits to both, compare the scores. Assign it to the
            #    category with the better score.
            elif smag_score > plant_score:
                f_smag.write(f"{read_id}\n")
                smag_count += 1
            elif plant_score > smag_score:
                f_plant.write(f"{read_id}\n")
                plant_count += 1
            # 4. If the scores are identical (a true tie), it's ambiguous.
            else:
                f_ambiguous.write(f"{read_id}\n")
                ambi_count += 1

    # Small summary to stdout (will land in PBS .o file)
    print(f"[classify_reads] alignments_seen={seen_alignments} "
          f"skipped_unmapped={skipped_unmapped} skipped_mapq={skipped_mapq} "
          f"skipped_flags={skipped_flags} plant={plant_count} smag={smag_count} ambiguous={ambi_count}")


# -------------------------------------------------------------------------
# PART 4: Script execution
# -------------------------------------------------------------------------
# This block checks if the script is being run directly from the command line.
if __name__ == "__main__":
    # Ensure the correct number of command-line arguments are provided.
    if len(sys.argv) != 4:
        print("Usage: python classify_reads.py <input_bam> <ref_list_dir> <output_dir>")
        sys.exit(1)
    
    # Assign command-line arguments to variables for clarity.
    input_bam_path = sys.argv[1]
    ref_list_dir = sys.argv[2]
    output_dir = sys.argv[3]
    
    # Call the main function to run the pipeline.
    main(input_bam_path, ref_list_dir, output_dir)