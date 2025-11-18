# Assembly Scripts

Scripts for metagenome assembly using metaMDBG and quality assessment.

## Prerequisites

- metaMDBG (v1.2+)
- minimap2 (v2.30+)
- MetaQUAST (QUAST v5.0+)
- MetaGeneMark
- R (v4.0+) with tidyverse, plotly
- Python 3.7+

## Workflow Overview
```
Filtered FASTQ
    ↓
[1] run_nanomdbg.sh → Raw contigs
    ↓
[2] polish_nanomdbg.sh → Polished contigs
    ↓
[3] assembly_qc.sh → Quality reports
    ↓
[4] extract_contig_headers.py → Contig metadata CSV
    ↓
[5] analyze_contigs.R → Statistics & filtered lists
```

---

## Configuration

Replace placeholders in scripts:

- `<PROJECT>`: NCI project code
- `<CONDA_ASSEMBLY_ENV>`: Path to conda environment
- `<METAMDBG_PATH>`: metaMDBG installation directory
- `<MINIMAP2_PATH>`: minimap2 binary directory
- `<METAGENEMARK_PATH>`: MetaGeneMark installation directory
- `<INPUT_FILTERED_FASTQ>`: Quality-filtered reads from step 02
- `<OUTPUT_ASSEMBLY_DIR>`: Directory for assembly output
- `<OUTPUT_POLISH_DIR>`: Directory for polished contigs
- `<OUTPUT_METAQUAST_DIR>`: Directory for QC reports

---

## Step 1: Assembly with metaMDBG

Assemble filtered reads using metaMDBG (nanoMDBG method for R10.4+ chemistry).
```bash
# Edit script
nano run_nanomdbg.sh

# Key parameters:
# --min-read-overlap 300  # Adjust based on read length (default: 1000)
# --min-abundance 1       # Include low-abundance k-mers

# Submit job
qsub run_nanomdbg.sh
```

**Output**: `<OUTPUT_ASSEMBLY_DIR>/contigs.fasta.gz`

**Resource requirements**:
- Queue: hugemem (required for large metagenomes)
- Memory: 1470 GB
- Walltime: 48 hours (may complete faster)
- Uses jobfs for faster I/O

---

## Step 2: Polish Contigs

Polish assembled contigs using original reads to improve accuracy.
```bash
# Edit script
nano polish_nanomdbg.sh

# Set paths:
# INPUT_CONTIGS: Output from step 1
# INPUT_READS: Original filtered reads (same as assembly input)

# Submit job
qsub polish_nanomdbg.sh
```

**Output**: `<OUTPUT_POLISH_DIR>/contigs.fasta.gz`

**Note**: Polishing typically improves base-level accuracy but is computationally expensive.

---

## Step 3: Quality Assessment

Assess assembly quality using MetaQUAST.
```bash
# Edit script
nano assembly_qc.sh

# Set paths:
# ASSEMBLY_FILES: Polished contigs from step 2
# NANOPORE_READS: Original reads

# Submit job
qsub assembly_qc.sh
```

**Output**: `<OUTPUT_METAQUAST_DIR>/metaquast_*/`
- `report.html`: Interactive quality report
- `report.txt`: Text summary
- Statistics on N50, misassemblies, genes, etc.

---

## Step 4: Extract Contig Metadata

Extract contig headers to CSV for downstream analysis.
```bash
# Run directly (fast, no PBS needed)
python extract_contig_headers.py \
    <INPUT_CONTIGS_FASTA> \
    <OUTPUT_HEADERS_CSV>

# Example:
python extract_contig_headers.py \
    /path/to/contigs.fasta.gz \
    /path/to/contig_headers.csv
```

**Output**: CSV with columns: `contig_id`, `length`, `coverage`, `circular`

---

## Step 5: Analyze and Filter Contigs

Generate statistics and filter contigs by quality criteria.
```bash
# Edit script
nano analyze_contigs.R

# Set paths:
# stats_file_path: CSV from step 4
# out_dir: Output directory for results

# Run
Rscript analyze_contigs.R
```

**Outputs**:
- `postassembly_qc_counts.csv`: Summary statistics
- `contigs_list_1Mb.txt`: Contigs ≥1 Mb
- `contigs_list_10kb_or_10x.txt`: Contigs ≥10 kb OR ≥10× coverage
- `contigs_list_20kb_or_20x.txt`: Contigs ≥20 kb OR ≥20× coverage
- `circular_contigs_list*.txt`: Circular contigs with various filters

**Filtering criteria**:
- High quality: ≥20 kb OR ≥20× coverage
- Medium quality: ≥10 kb OR ≥10× coverage
- Large contigs: ≥1 Mb
- Circular: Likely complete genomes/plasmids

---

## Expected Results

Typical metaMDBG assembly metrics:
- **N50**: 50-200 kb (depends on community complexity)
- **Total assembly size**: Varies widely (100 Mb - 10 Gb)
- **Circular contigs**: 10-100 (putative complete genomes)
- **Longest contig**: Often >1 Mb for dominant species

**Good assembly indicators**:
- High N50 (>100 kb)
- Many circular contigs
- Low misassembly rate in MetaQUAST
- High gene completeness

---

## Troubleshooting

### Assembly Too Fragmented

**Symptoms**: Low N50, many short contigs

**Solutions**:
1. Increase read quality threshold in step 01-quality-control
2. Adjust `--min-read-overlap` (try 500 instead of 300)
3. Check read length distribution (may need longer reads)

---

### Out of Memory

**Symptoms**: Job killed, "out of memory" errors

**Solutions**:
1. Assembly already uses hugemem queue (max available)
2. Reduce dataset size (subsample reads)
3. Increase `--min-abundance` to filter rare k-mers

---

### Very Long Runtime

**Symptoms**: Job approaching 48-hour limit

**Solutions**:
1. Consider using `--max-k` to limit iterations
2. Increase `--min-abundance` to speed up
3. Split large datasets and assemble separately

---

### No Circular Contigs

**Symptoms**: `circular=no` for all contigs

**Possible causes**:
- Low coverage
- High diversity (strain variation)
- Incomplete sequencing

**Not necessarily bad**: Many metagenomes don't have complete genomes

---

## Parameter Tuning Guide

### For Short Reads (N50 < 5 kb)
```bash
--min-read-overlap 300   # Lower overlap requirement
--min-abundance 2        # Filter more noise
```

### For Long Reads (N50 > 10 kb)
```bash
--min-read-overlap 1000  # Default, more stringent
--min-abundance 1        # Keep rare sequences
```

### For High Diversity Communities
```bash
--min-abundance 2        # Reduce complexity
--max-k 10               # Limit iterations
```

### For Low Diversity (e.g., enrichments)
```bash
--min-abundance 1        # Keep all sequences
--max-k 0                # No limit (default)
```

---

## References

- metaMDBG: https://github.com/GaetanBenoitDev/metaMDBG
- MetaQUAST: http://quast.sourceforge.net/metaquast
- MetaGeneMark: http://exon.gatech.edu/meta_gmhmmp.cgi