#!/bin/bash
#PBS -P <PROJECT>
#PBS -q hugemem
#PBS -l walltime=15:00:00
#PBS -l mem=1000GB
#PBS -l jobfs=1000GB
#PBS -l ncpus=48
#PBS -l storage=scratch/<PROJECT>+gdata/<PROJECT>
#PBS -l wd

set -euo pipefail
set -x

# ============================================================================
# Assembly QC - MetaQUAST Quality Assessment
# ============================================================================
# Description: Assess assembly quality using MetaQUAST
# Input: Assembled contigs (FASTA) and Nanopore reads (FASTQ)
# Output: MetaQUAST reports and statistics
# ============================================================================

echo "Job started on: $(date)"
echo "Host: $(hostname)"

# Activate conda environment
source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate <CONDA_ASSEMBLY_ENV>

# Add MetaGeneMark to PATH
export PATH="<METAGENEMARK_PATH>/mgm:$PATH"

# Input files
ASSEMBLY_FILES=(
    "<INPUT_ASSEMBLY_FASTA>"
)

NANOPORE_READS="<INPUT_NANOPORE_READS>"
BASE_OUTPUT_DIR="<OUTPUT_METAQUAST_DIR>"

# Process each assembly
for ASSEMBLY in "${ASSEMBLY_FILES[@]}"; do
    # Extract assembler and dataset names
    ASSEMBLER_DIR_NAME=$(basename "$(dirname "$(dirname "$ASSEMBLY")")")
    ASSEMBLER=$(echo "$ASSEMBLER_DIR_NAME" | cut -d'_' -f1)
    
    DATASET=$(basename "$(dirname "$ASSEMBLY")" | sed 's/_FASTQ_filtered.*//;s/_myloasm_assembly//')
    OUTPUT_DIR="${BASE_OUTPUT_DIR}/metaquast_${DATASET}_${ASSEMBLER}/"
    mkdir -p "$OUTPUT_DIR"
    
    ASSEMBLY_BASENAME=$(basename "$ASSEMBLY")
    ASSEMBLY_LABEL="${DATASET}_${ASSEMBLER}"
    
    # Copy files to jobfs for faster processing
    LOCAL_ASSEMBLY="$PBS_JOBFS/$ASSEMBLY_BASENAME"
    LOCAL_READS="$PBS_JOBFS/$(basename "$NANOPORE_READS")"
    LOCAL_OUTPUT="$PBS_JOBFS/metaquast_out_${ASSEMBLY_LABEL}"
    mkdir -p "$LOCAL_OUTPUT"
    
    echo "Copying files to jobfs..."
    cp "$ASSEMBLY" "$LOCAL_ASSEMBLY"
    cp "$NANOPORE_READS" "$LOCAL_READS"
    
    echo "Running MetaQUAST on: $ASSEMBLY_LABEL"
    
    # MetaQUAST quality assessment
    #
    # Basic options:
    #   --nanopore              Nanopore reads file
    #   --labels                Assembly label
    #   --output-dir            Output directory
    #   --threads               Number of threads
    #
    # Alignment options:
    #   --min-contig            Minimum contig length [1000]
    #   --min-identity          Minimum identity for alignment [90.0]
    #   --min-alignment         Minimum alignment length [200]
    #   --ambiguity-usage       Use all alignments [all]
    #   --ambiguity-score       Ambiguity score threshold [0.95]
    #
    # Misassembly detection:
    #   --extensive-mis-size    Extensive misassembly size [1000]
    #   --scaffold-gap-max-size Max scaffold gap size [10000]
    #   --unaligned-part-size   Unaligned part size [200]
    #
    # Feature detection:
    #   --x-for-Nx              Calculate NX for X% [90]
    #   --rna-finding           Enable RNA gene finding
    #   --gene-finding          Enable gene finding
    #   --conserved-genes-finding Enable conserved gene finding
    #   --mgm                   Use MetaGeneMark for gene finding
    #   --max-ref-number        Max reference genomes [0 = no references]
    
    <CONDA_ASSEMBLY_ENV>/bin/metaquast.py \
        "$LOCAL_ASSEMBLY" \
        --nanopore "$LOCAL_READS" \
        --labels "$ASSEMBLY_LABEL" \
        --min-contig 1000 \
        --min-identity 90.0 \
        --min-alignment 200 \
        --ambiguity-usage all \
        --ambiguity-score 0.95 \
        --extensive-mis-size 1000 \
        --scaffold-gap-max-size 10000 \
        --unaligned-part-size 200 \
        --x-for-Nx 90 \
        --rna-finding \
        --gene-finding \
        --conserved-genes-finding \
        --mgm \
        --max-ref-number 0 \
        --output-dir "$LOCAL_OUTPUT" \
        --threads "$PBS_NCPUS"
    
    # Copy results to permanent storage
    echo "Copying results to permanent storage..."
    cp -r "$LOCAL_OUTPUT"/* "$OUTPUT_DIR"/
done

echo "MetaQUAST processing complete for all assemblies."
echo "Job completed on: $(date)"