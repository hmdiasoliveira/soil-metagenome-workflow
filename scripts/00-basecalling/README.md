# Basecalling Scripts

PBS scripts for Oxford Nanopore basecalling using Dorado on GPU nodes.

## Prerequisites

- Dorado basecaller (v1.1.0 or later)
- CUDA-capable GPU (V100 or A100)
- Access to NCI GPU queues (`gpuvolta` or `dgxa100`)
- Modules: `cuda`, `samtools`, `parallel`

## Workflow Overview
```
POD5 files
    ↓
[1] run_dorado_basecaller.sh → Basecalled BAM (all barcodes)
    ↓
[2] run_dorado_demux.sh → Per-barcode BAM + FASTQ.gz files
```

If basecalling fails partway through:
```
[1b] resume_dorado_basecaller.sh → Resume from incomplete BAM
```

---

## Configuration

Before running, edit each script to replace placeholders:

### Common Placeholders

- `<PROJECT>`: Your NCI project code (e.g., `xe2`)
- `<DORADO_INSTALLATION_PATH>`: Path to Dorado installation directory
  - Example: `/g/data/xe2/adria/apps/dorado-1.1.0-linux-x64`
- `<DORADO_MODEL_PATH>`: Path to Dorado model directory
  - Example: `/g/data/xe2/adria/apps/dorado_persistent_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0`
- `<BARCODE_KIT>`: Barcoding kit used for library preparation
  - Example: `SQK-NBD114-96` (96-sample barcoding kit)

### Script-Specific Placeholders

**run_dorado_basecaller.sh**:
- `<INPUT_POD5_DIR>`: Directory containing POD5 files
- `<OUTPUT_BAM_DIR>`: Directory for output BAM file

**resume_dorado_basecaller.sh**:
- `<INPUT_POD5_DIR>`: Directory containing POD5 files (same as original run)
- `<OUTPUT_BAM_DIR>`: Directory for output BAM file
- `<INCOMPLETE_BAM_PATH>`: Path to incomplete BAM from failed run
- `<RESUMED_BAM_OUTPUT_PATH>`: Path for new, completed BAM file

**run_dorado_demux.sh**:
- `<INPUT_BAM_PATH>`: Path to basecalled BAM file (or directory)
- `<OUTPUT_DEMUX_BAM_DIR>`: Directory for per-barcode BAM files
- `<OUTPUT_FASTQ_DIR>`: Directory for per-barcode FASTQ files

---

## Step 1: Initial Basecalling

Basecall POD5 files with inline demultiplexing and adapter trimming.
```bash
# Edit script to set paths
nano run_dorado_basecaller.sh

# Key settings:
# - GPU queue: gpuvolta (V100) or dgxa100 (A100)
# - Model: dna_r10.4.1_e8.2_400bps_sup@v5.2.0 (super-accurate)
# - Barcode kit: SQK-NBD114-96 (or your kit)

# Submit job
qsub run_dorado_basecaller.sh
```

**Output**: `<OUTPUT_BAM_DIR>/<input_dir_name>_trimmed.bam`

**Options explained**:
- `--device cuda:all`: Use all available GPUs
- `--kit-name`: Barcode kit for demultiplexing
- `--trim all`: Remove adapters and barcodes from reads
- `--recursive`: Process POD5 files in subdirectories

---

## Step 2: Resume Incomplete Basecalling (If Needed)

If basecalling job timed out or failed, resume from the incomplete BAM.
```bash
# Edit script to set paths
nano resume_dorado_basecaller.sh

# Set paths:
# - incompleteBam: Output from failed run
# - resumedBam: New output file name
# - Same InputDir, modelCfg, barcodeKit as original

# Submit job
qsub resume_dorado_basecaller.sh
```

**Output**: `<RESUMED_BAM_OUTPUT_PATH>`

**Notes**:
- Uses same POD5 input directory as original run
- Dorado automatically skips already-processed reads
- Can resume multiple times if needed

---

## Step 3: Demultiplexing and FASTQ Conversion

Separate barcodes into individual files and convert BAM to FASTQ.
```bash
# Edit script to set paths
nano run_dorado_demux.sh

# Set paths:
# - inputBAM: Output from basecalling step
# - demuxDir: Per-barcode BAM files
# - fastqDir: Per-barcode FASTQ files (gzipped)

# Submit job
qsub run_dorado_demux.sh
```

**Output**: 
- `<OUTPUT_DEMUX_BAM_DIR>/barcode*.bam` (one per barcode)
- `<OUTPUT_FASTQ_DIR>/barcode*.fastq.gz` (one per barcode)

**Options explained**:
- `--no-classify`: Skip ML classification (use kit-based demux)
- `--emit-fastq`: Also output FASTQ files directly
- `parallel -j 12`: Convert 12 BAM files to FASTQ simultaneously

---

## Resource Requirements

### run_dorado_basecaller.sh (V100 GPU)

| Resource    | Value     | Notes                                |
|-------------|-----------|--------------------------------------|
| Queue       | gpuvolta  | V100 GPU nodes                       |
| GPUs        | 1         | Single GPU sufficient                |
| CPUs        | 12        | For data I/O                         |
| Memory      | 384 GB    | Large for complex libraries          |
| Walltime    | 48 hours  | Adjust based on data size            |
| Jobfs       | 400 GB    | Temporary storage                    |

**Typical performance**: ~500 kb/s per GPU (V100)

---

### resume_dorado_basecaller.sh (A100 GPU)

| Resource    | Value     | Notes                                |
|-------------|-----------|--------------------------------------|
| Queue       | dgxa100   | A100 GPU nodes (faster)              |
| GPUs        | 1         | Single GPU                           |
| CPUs        | 16        | Higher for faster I/O                |
| Memory      | 384 GB    | Same as initial run                  |
| Walltime    | 48 hours  | Conservative estimate                |
| Jobfs       | 1200 GB   | Larger for resume buffer             |

**Typical performance**: ~1-2 Mb/s per GPU (A100)

---

### run_dorado_demux.sh

| Resource    | Value     | Notes                                |
|-------------|-----------|--------------------------------------|
| Queue       | dgxa100   | Can use GPU queue but doesn't need GPU for demux |
| GPUs        | 1         | Not used, but required for queue     |
| CPUs        | 16        | Parallel BAM→FASTQ conversion        |
| Memory      | 96 GB     | Moderate for demux                   |
| Walltime    | 20 hours  | Usually completes faster             |
| Jobfs       | 100 GB    | Minimal temporary storage            |

**Notes**: 
- Demux step is CPU-bound, not GPU-bound
- Can alternatively run on normal queue without GPU

---

## Barcode Kit Reference

Common Oxford Nanopore barcoding kits:

| Kit Name         | Barcodes | Notes                               |
|------------------|----------|-------------------------------------|
| SQK-NBD114-96    | 1-96     | 96-sample multiplexing              |
| SQK-NBD114-24    | 1-24     | 24-sample multiplexing              |
| SQK-RBK114-96    | 1-96     | Rapid barcoding (96 samples)        |
| SQK-RBK114-24    | 1-24     | Rapid barcoding (24 samples)        |
| SQK-16S114-24    | 1-24     | 16S amplicon barcoding              |

Specify exactly as shown in Dorado's `--kit-name` option.

---

## Model Selection

### Recommended Models (R10.4.1 chemistry)

| Model                                      | Accuracy | Speed     | Use Case              |
|--------------------------------------------|----------|-----------|-----------------------|
| dna_r10.4.1_e8.2_400bps_sup@v5.2.0        | Highest  | Slow      | Publication-quality   |
| dna_r10.4.1_e8.2_400bps_hac@v5.2.0        | High     | Medium    | Standard analysis     |
| dna_r10.4.1_e8.2_400bps_fast@v5.2.0       | Moderate | Fast      | Quick QC              |

**sup** (super-accurate) is recommended for metagenomics to minimize assembly errors.

---

## Troubleshooting

### Job Killed or Timed Out

**Symptoms**: Job stops before completion, incomplete BAM file

**Solution**: Use `resume_dorado_basecaller.sh`
1. Note the incomplete BAM file path
2. Edit resume script with paths
3. Submit resume job
4. Repeat if necessary

### Nested Directories from Demux

After running `run_dorado_demux.sh`, you get nested directories instead of per-barcode files:
```
demux_BAM/
├── Batch3_I/
│   └── 20251029_0909_0_PBE46605_ba040519/
│       └── fastq_pass/
│           └── PBE46605_pass_ba040519_00000000_0.fastq
```

Instead of the expected:
```
demux_BAM/
├── barcode01.fastq
├── barcode02.fastq
...
```

### Cause

This happens when:
1. You merged multiple BAM files from different sequencing runs
2. The BAM has complex read group (RG) structure
3. Dorado demux organizes by run metadata instead of just barcodes

### Solution: Use `split_fastq_by_barcode.sh`

Use this script when:
1. You already have a mixed FASTQ file from dorado demux
2. The FASTQ contains reads from multiple barcodes mixed together
3. Barcode information is in the read headers (RG tags)
4. You want to separate them without going back to BAM

### Usage
```bash
# Edit script to set paths
nano split_fastq_by_barcode.sh

# Set:
# - inputFASTQ: Path to mixed FASTQ file
# - outputDir: Output directory for per-barcode files

# Submit job
qsub split_fastq_by_barcode.sh
```

### Example

**Input FASTQ header:**
```
@c34bef41-4470-4c94-9510-b4ddbd066df6 RG:Z:ba040519...SQK-NBD114-96_barcode11
```

**Output:**
```
outputDir/
├── barcode01.fastq.gz
├── barcode11.fastq.gz
├── barcode05.fastq.gz
...
└── unclassified.fastq.gz
```

### Features

- Handles both compressed (.gz) and uncompressed FASTQ
- Zero-pads barcode numbers (barcode01, barcode02, etc.)
- Reports read counts per barcode
- Automatically creates output files as needed
- Memory efficient (streams data)

---

### CUDA Out of Memory

**Symptoms**: `CUDA error: out of memory`

**Solutions**:
1. Reduce `--chunk-size` (not in current scripts, but can add)
2. Use fewer GPUs: `--device cuda:0` instead of `cuda:all`
3. Request more GPU memory (A100 has more than V100)

---

### No Reads Demultiplexed

**Symptoms**: Empty barcode BAM files after demux

**Possible causes**:
- Incorrect barcode kit specified
- Reads not barcoded (check library prep)
- Barcodes removed by `--trim all` but not demuxed

**Solution**: 
- Verify barcode kit matches library prep
- Check basecalling log for demux statistics
- Re-run basecalling if wrong kit was used

---

### Slow Basecalling Speed

**Symptoms**: Basecalling taking much longer than expected

**Optimization**:
1. Use A100 GPUs (`dgxa100` queue) instead of V100
2. Ensure model is pre-downloaded (set `DORADO_MODELS_DIRECTORY`)
3. Check POD5 file integrity
4. Increase `--chunk-size` if enough GPU memory

---

## Output File Structure

After completing all steps:
```
<OUTPUT_BAM_DIR>/
└── <input_name>_trimmed.bam              # All barcodes, basecalled

<OUTPUT_DEMUX_BAM_DIR>/
├── unclassified.bam                       # Reads without barcode
├── barcode01.bam
├── barcode02.bam
├── ...
└── barcode96.bam

<OUTPUT_FASTQ_DIR>/
├── unclassified.fastq.gz
├── barcode01.fastq.gz
├── barcode02.fastq.gz
├── ...
└── barcode96.fastq.gz
```

---

## Expected Output Statistics

After demultiplexing, check barcode distribution:
```bash
# Count reads per barcode
for f in <OUTPUT_FASTQ_DIR>/*.fastq.gz; do
    echo "$(basename $f): $(zcat $f | grep -c '^@')"
done
```

**Typical distribution**:
- Most barcodes: 1-5M reads each
- Unclassified: <5% of total reads
- Empty barcodes: Indicate unused barcodes or failed samples

---

## Integration with QC Workflow

After basecalling and demultiplexing, proceed to quality control:
```bash
# 1. Complete basecalling workflow
qsub scripts/00-basecalling/run_dorado_basecaller.sh
qsub scripts/00-basecalling/run_dorado_demux.sh

# 2. Run QC on demultiplexed FASTQ files
# Edit INPUT_DIR to point to <OUTPUT_FASTQ_DIR>
qsub scripts/01-quality-control/run_nanoplot.sh
qsub scripts/01-quality-control/run_fastqc.sh
qsub scripts/01-quality-control/multiqc_report.sh

# 3. Filter based on QC results
qsub scripts/01-quality-control/filter_chopper.sh
```

---

## Example: Complete Basecalling Workflow
```bash
# Step 1: Initial basecalling (48 hours on V100)
qsub run_dorado_basecaller.sh

# Check job status
qstat -u $USER

# If job completes successfully:
ls <OUTPUT_BAM_DIR>/*.bam

# Step 2: Demultiplex (4-8 hours)
qsub run_dorado_demux.sh

# Verify outputs
ls <OUTPUT_FASTQ_DIR>/*.fastq.gz | wc -l  # Should match number of barcodes used

# Step 3: Proceed to quality control
cd ../01-quality-control
```

---

## Advanced: Inline vs. Post-Demultiplexing

### Inline Demultiplexing (Current Workflow)
```bash
dorado basecaller MODEL INPUT --kit-name KIT --trim all > output.bam
dorado demux --output-dir DEMUX_DIR output.bam
```

**Advantages**:
- Single basecalling pass
- Faster overall
- Easier to resume if interrupted

**Disadvantages**:
- All barcodes in one large file initially
- Harder to process individual samples early

---

### Post-Demultiplexing (Alternative)
```bash
# Basecall without demux
dorado basecaller MODEL INPUT --trim all > output.bam

# Demux later with different kit if needed
dorado demux --kit-name KIT --output-dir DEMUX_DIR output.bam
```

**Use case**: When barcode kit is uncertain or might change

---

## Dorado Version Notes

These scripts were tested with **Dorado v1.1.0**.

For other versions:
- Check `dorado basecaller --help` for option changes
- Model paths may differ
- Resume functionality may not be available in older versions

---

## GPU Queue Access

NCI GPU queues require special project access:

- **gpuvolta**: V100 GPUs (older, more available)
- **dgxa100**: A100 GPUs (faster, more in-demand)

If you don't have access:
```bash
# Check project GPU quota
nci_account -P <PROJECT>

# Request GPU access via NCI helpdesk if needed
```

---

## References

- Dorado documentation: https://github.com/nanoporetech/dorado
- Oxford Nanopore community: https://community.nanoporetech.com/
- NCI GPU documentation: https://opus.nci.org.au/display/Help/GPU+Nodes