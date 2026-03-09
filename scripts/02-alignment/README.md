# Contamination Removal Scripts

PBS scripts for identifying and removing plant host contamination from soil metagenomic data using competitive mapping.

## Overview

This workflow uses **competitive mapping** to distinguish microbial reads from plant host contamination:

1. **Competitive Mapping**: Map reads to combined plant+SMAG reference
2. **Read Classification**: Compare alignment scores to categorize each read
3. **Contamination Filtering**: Extract clean microbial FASTQ files

## Prerequisites

- Minimap2 (v2.30+)
- Samtools (v1.19+)
- Python 3.7+ with pysam
- seqtk (for FASTQ filtering)
- Combined reference genome (plant + SMAG genomes)
- Reference contig lists (plant_contigs.txt, smag_contigs.txt)

## Workflow Diagram
```
Filtered FASTQ
    ↓
[1] competitive_mapping.sh → Sorted BAM files
    ↓
[2] classify_reads.py → Read name lists (plant/SMAG/ambiguous)
    ↓
[3] filter_contaminants.sh → Clean microbial FASTQ files
```

---

## Preserving BAM Tags (MM/ML) Through FASTQ Conversion

When working with modified basecalled data (e.g. Dorado with `--modified-bases 5mC_5hmC`), unaligned BAMs contain MM and ML tags encoding per-read methylation probabilities. These tags are **lost** during standard BAM → FASTQ conversion because FASTQ format does not natively carry SAM auxiliary tags.

### Why This Matters

If your input BAMs contain modification calls (MM/ML tags) and you need to preserve them through alignment, a standard `samtools fastq` → `minimap2` → `samtools sort` pipeline will silently drop the methylation information.

### Option 1: `samtools fastq -T` + `minimap2 -y`

Use `samtools fastq -T` to carry specified tags into the FASTQ header comments, then `minimap2 -y` to copy those comments into the SAM output:

```bash
samtools fastq -T MM,ML unaligned.bam \
  | minimap2 -y -ax map-ont reference.fa - \
  | samtools sort -o aligned.bam
```

- `-T MM,ML`: writes the MM and ML tags into the FASTQ header comment field
- `-y`: copies FASTQ header comments into the SAM auxiliary tags

To carry **all** tags (not just MM/ML):
```bash
samtools fastq -T '*' unaligned.bam \
  | minimap2 -y -ax map-ont reference.fa - \
  | samtools sort -o aligned.bam
```

### Option 2: `dorado aligner`

Dorado's built-in aligner accepts unaligned BAM directly (any HTS format) and preserves MM/ML tags natively, bypassing FASTQ conversion entirely:

```bash
dorado aligner reference.fa unaligned.bam > aligned.bam
```

Notes on `dorado aligner`:
- Uses the `lr:hq` minimap2 preset by default (suited for sup-basecalled high-accuracy reads)
- Override with `--mm2-opts "-x map-ont"` if needed
- Supports `--output-dir` for batch processing

### Verifying Tags Are Preserved

After alignment, confirm MM/ML tags are present in the output:
```bash
samtools view aligned.bam | head -1 | tr '\t' '\n' | grep "^MM\|^ML"
```

If no output is returned, the tags were lost during conversion/alignment.

### Relevance to This Workflow

The competitive mapping step in this pipeline converts FASTQ to BAM via minimap2. If your upstream FASTQ files were generated from modified basecalled BAMs **without** `-T MM,ML`, the methylation tags will already be absent. To preserve them:

1. Start from the original unaligned BAMs (not FASTQ)
2. Use `samtools fastq -T MM,ML` piped into `minimap2 -y` for the competitive mapping step
3. Or use `dorado aligner` for the initial alignment and then filter the resulting BAM

---

## Configuration

### Common Placeholders

- `<PROJECT>`: Your NCI project code
- `<MINIMAP2_PATH>`: Path to minimap2 binary directory
- `<CONDA_PATH>`: Path to conda installation
- `<ALIGN_ENV_PATH>`: Path to alignment conda environment
- `<SCRIPT_DIR>`: Directory containing Python scripts

### Reference Files

**Combined Reference**:
- `<COMBINED_REFERENCE_MMI>`: Minimap2 index (.mmi file)
- `<COMBINED_REFERENCE_DICT>`: Sequence dictionary (.dict file)

**Contig Lists**:
- `<REFERENCE_CONTIG_LISTS_DIR>`: Directory containing:
  - `plant_contigs.txt`: One plant contig name per line
  - `smag_contigs.txt`: One SMAG contig name per line

---

## Step 1: Build Combined Reference

Before running the workflow, prepare the combined reference:
```bash
# Concatenate plant and SMAG genomes
cat plant_genomes/*.fasta smag_catalogue.fasta > combined_reference.fasta

# Build minimap2 index
minimap2 -d combined_reference.mmi combined_reference.fasta

# Create sequence dictionary for samtools
samtools dict combined_reference.fasta > combined_reference.dict

# Extract contig names
grep "^>" plant_genomes/*.fasta | sed 's/>//' | cut -f1 -d' ' > plant_contigs.txt
grep "^>" smag_catalogue.fasta | sed 's/>//' | cut -f1 -d' ' > smag_contigs.txt
```

---

## Step 2: Competitive Mapping

Map reads to the combined reference to enable competitive alignment scoring.
```bash
# Edit script to set paths
nano competitive_mapping.sh

# Key settings:
# - REFERENCE_MMI: Path to combined reference index
# - DICT_FILE: Path to sequence dictionary
# - INPUT_DIR: Directory with filtered FASTQ files
# - OUTPUT_DIR: Output directory for BAM files

# Submit job
qsub competitive_mapping.sh
```

**Output**: `<OUTPUT_DIR>/<sample>.sorted.bam` (plus .bai and .flagstat.txt)

**Mapping Options**:
- `-ax map-ont`: Oxford Nanopore preset
- `--MD`: Include MD tag for variant calling
- `--secondary=yes`: Retain secondary alignments (default)
- `-N 5`: Keep up to 5 secondary alignments
- `-p 0.8`: Secondary must have ≥80% of primary score

---

## Step 3: Classify Reads

Use Python script to classify reads based on best alignment scores.

### 3a. Run Classification Independently (Optional)
```bash
# Activate conda environment with pysam
conda activate <ALIGN_ENV_PATH>

# Run for a single sample
python classify_reads.py \
    input.sorted.bam \
    /path/to/reference_contig_lists \
    output_dir

# Outputs:
#   output_dir/plant_unique_names.txt
#   output_dir/smag_unique_names.txt
#   output_dir/ambiguous_names.txt
```

### 3b. Or Use Integrated Workflow (Recommended)

Skip to Step 4 - classification is integrated into `filter_contaminants.sh`.

---

## Step 4: Filter Contaminants (Integrated Workflow)

Complete pipeline: classify reads, split BAMs, calculate statistics, and create clean FASTQ.
```bash
# Edit script to set paths
nano filter_contaminants.sh

# Required paths:
# - BAM_DIR: Output from competitive_mapping.sh
# - REF_LIST_DIR: Directory with plant_contigs.txt and smag_contigs.txt
# - FASTQ_DIR: Original filtered FASTQ files
# - SCRIPT_DIR: Directory containing classify_reads.py

# Submit job
qsub filter_contaminants.sh
```

**This script performs**:
1. Read classification (calls classify_reads.py)
2. BAM splitting (plant/SMAG/ambiguous)
3. Contamination statistics calculation
4. Clean FASTQ generation (SMAG + unmapped reads)

**Outputs**:
```
<BAM_DIR>/
├── split_reads/
│   └── <sample>/
│       ├── plant_unique_names.txt
│       ├── smag_unique_names.txt
│       ├── ambiguous_names.txt
│       ├── plant_unique_final.bam
│       ├── smag_unique_final.bam
│       └── ambiguous_final.bam
├── idxstats/
│   └── <sample>.idxstats.txt
└── final_output/
    ├── <sample>_smag_and_unmapped_clean.fastq.gz
    └── all_barcodes_read_proportions.txt
```

---

## Classification Logic

### Scoring System

For each read, the script compares alignments using a tuple:
```
(alignment_score, MAPQ, aligned_length)
```

Tuples are compared element-by-element (AS first, then MAPQ, then length).

### Decision Rules

| Condition | Category | Rationale |
|-----------|----------|-----------|
| Only plant hits | Plant-unique | No microbial similarity |
| Only SMAG hits | SMAG-unique | No plant similarity |
| Both, SMAG score > plant score | SMAG-unique | Better microbial match |
| Both, plant score > SMAG score | Plant-unique | Better plant match |
| Both, scores exactly equal | Ambiguous | True tie |

### Handling Edge Cases

**Unmapped reads**: 
- Not analyzed by classify_reads.py
- Added to clean FASTQ in filter_contaminants.sh
- Rationale: Likely novel microbes not in reference

**Secondary alignments**:
- Included by default (`INCLUDE_SECONDARY = True`)
- Important for ONT reads that span multiple genes
- Ensures best hit is found even if not primary

**Missing alignment scores**:
- AS defaults to -10⁹ if absent
- Ensures real scores always win ties

---

## Resource Requirements

### competitive_mapping.sh

| Resource | Value    | Notes                              |
|----------|----------|------------------------------------|
| Queue    | hugemem  | Large reference requires RAM       |
| Memory   | 1000 GB  | Combined plant+SMAG reference      |
| CPUs     | 48       | Parallel mapping                   |
| Walltime | 48 hours | Conservative for large datasets    |
| Jobfs    | 500 GB   | Copy reference and inputs locally  |

**Performance**: ~5-10 minutes per Gbp (depends on reference size)

---

### filter_contaminants.sh

| Resource | Value    | Notes                              |
|----------|----------|------------------------------------|
| Queue    | hugemem  | Multiple concurrent operations     |
| Memory   | 1000 GB  | Handles large BAM files            |
| CPUs     | 48       | Parallel BAM operations            |
| Walltime | 48 hours | Full pipeline for many samples     |
| Jobfs    | 500 GB   | Temporary BAM files                |

**Note**: Uses half available CPUs per samtools operation to avoid oversubscription.

---

## Understanding Output Statistics

### Contamination Report

`all_barcodes_read_proportions.txt` shows percentage of reads mapping to each reference:
```
Barcode    Total_Mapped_Reads  plant_reads(%)  smag_reads(%)
barcode01  5234567             1234567(23.6%)  4000000(76.4%)
barcode02  6789012             567890(8.4%)    6221122(91.6%)
```

**Interpretation**:

| Plant % | Microbial % | Sample Type | Action |
|---------|-------------|-------------|--------|
| >50%    | <50%        | Root tissue | Poor sampling |
| 20-40%  | 60-80%      | Rhizosphere | Expected |
| 5-15%   | 85-95%      | Bulk soil   | Good |
| <5%     | >95%        | Pure culture| Excellent |

### Microbe:Plant Ratio

Calculate the enrichment of microbial reads:
```bash
# From proportions file
awk -F'\t' 'NR>1 {
    split($3, plant, "[(:]")
    split($4, smag, "[(:]")
    ratio = smag[2] / plant[2]
    print $1, ratio
}' all_barcodes_read_proportions.txt
```

**Expected ratios**:
- Rhizosphere: 3:1 to 10:1
- Bulk soil: >20:1

---

## Troubleshooting

### Issue 1: SAM File Missing @SQ Headers

**Symptoms**: Samtools fails with "missing @SQ header"

**Cause**: Some minimap2 versions omit @SQ lines when using prebuilt .mmi

**Solution**: Script automatically adds headers from .dict file
```bash
if ! grep -q "^@SQ" "$SAM_FILE"; then
    cat "$DICT_FILE" "$SAM_FILE" > "${SAM_FILE}_fixed.sam"
fi
```

**To create .dict manually**:
```bash
samtools dict reference.fasta > reference.dict
```

---

### Issue 2: High Memory Usage During Filtering

**Symptoms**: Job killed by OOM (out of memory)

**Solutions**:
1. Reduce concurrent operations (use fewer threads)
2. Process samples sequentially instead of in parallel
3. Use two-pass filtering (already implemented):
```bash
   # Pass 1: Filter by read names
   samtools view -N names.txt input.bam > temp.bam
   # Pass 2: Filter by region
   samtools view -L regions.bed temp.bam > output.bam
```

---

### Issue 3: Empty SMAG BAM Files

**Symptoms**: `smag_unique_final.bam` has no reads

**Possible causes**:
1. All reads are plant-derived (check proportions file)
2. Contig name mismatch between BAM and `smag_contigs.txt`
3. All reads classified as ambiguous

**Diagnosis**:
```bash
# Check read counts
wc -l split_reads/*/plant_unique_names.txt
wc -l split_reads/*/smag_unique_names.txt
wc -l split_reads/*/ambiguous_names.txt

# Check contig name format
head smag_contigs.txt
samtools view -H input.bam | grep "^@SQ" | head
```

**Fix**: Ensure contig names match exactly (no spaces, special characters consistent)

---

### Issue 4: Python Script Crashes with Memory Error

**Symptoms**: `classify_reads.py` killed during execution

**Cause**: Very large BAM files exhaust RAM while building dictionaries

**Solutions**:
1. Increase job memory (already at 1 TB)
2. Process BAM in chunks (requires script modification)
3. Use samtools to pre-filter unmapped reads:
```bash
   samtools view -F 4 -b input.bam > mapped_only.bam
   python classify_reads.py mapped_only.bam ...
```

---

### Issue 5: Slow seqtk Extraction

**Symptoms**: Task 3 takes many hours

**Cause**: Large read name lists + gzipped FASTQ

**Optimization**:
```bash
# Sort read names for faster seqtk lookup
sort smag_reads_to_keep.txt > smag_reads_sorted.txt

# Use parallel compression
seqtk subseq input.fastq.gz smag_reads_sorted.txt | pigz -p 8 > output.fastq.gz
```

---

## Advanced: Adjusting Classification Thresholds

Edit `classify_reads.py` to modify filtering behavior:

### Require Minimum MAPQ
```python
# Line ~17
MIN_MAPQ = 20  # Only consider alignments with MAPQ ≥ 20
```

**Use case**: Reduce ambiguous reads by requiring confident mappings

---

### Exclude Secondary Alignments
```python
# Line ~22
INCLUDE_SECONDARY = False  # Only use primary alignments
```

**Use case**: Simpler classification logic, faster execution

**Trade-off**: May miss best hit if it's secondary

---

### Custom Scoring Function
```python
# Line ~29
def _score_tuple(read):
    # Weight MAPQ more heavily
    ascore = read.get_tag('AS') if read.has_tag('AS') else -10**9
    mapq = (read.mapping_quality or 0) * 10  # 10× weight
    alen = read.query_alignment_length or 0
    return (ascore, mapq, alen)
```

**Use case**: Prioritize uniquely mapped reads (high MAPQ)

---

## Integration with Assembly Workflow

After contamination removal, proceed to assembly:
```bash
# 1. Complete contamination removal
qsub scripts/02-contamination-removal/competitive_mapping.sh
qsub scripts/02-contamination-removal/filter_contaminants.sh

# 2. Verify clean FASTQ files
ls final_output/*_smag_and_unmapped_clean.fastq.gz

# 3. Proceed to assembly
# Edit INPUT_DIR in assembly script to point to final_output/
qsub scripts/03-assembly/run_nanomdbg.sh
```

---

## Expected Timeline

For a typical 96-sample multiplexed run:

| Step | Duration | Notes |
|------|----------|-------|
| Competitive mapping | 24-36 hours | Depends on reference size |
| Classification | 4-8 hours | Python script (per sample) |
| Filtering | 12-24 hours | BAM operations + FASTQ extraction |
| **Total** | **2-3 days** | Sequential execution |

**Optimization**: Run mapping and filtering jobs in parallel for independent samples.

---

## Validation Checklist

After completing the workflow:

- [ ] All samples have BAM files in `split_reads/`
- [ ] Proportions file shows reasonable contamination levels
- [ ] Clean FASTQ files are not empty
- [ ] Read counts match: FASTQ ≈ (SMAG reads + unmapped reads)
- [ ] Ambiguous reads < 5% of total
- [ ] Plant-unique reads match expected contamination level

---

## References

- Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100.
- Samtools: Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008.
- Competitive mapping strategy: Adapted from metagenome host decontamination protocols in HMP and MGnify.