# Basecalling Scripts

PBS scripts for Oxford Nanopore basecalling using Dorado on GPU nodes.

## Prerequisites

- Dorado basecaller (v1.1.0 or later)
- CUDA-capable GPU (V100 or A100)
- Access to NCI GPU queues (`gpuvolta` or `dgxa100`)
- Modules: `cuda`, `samtools`, `parallel`

---

## Workflow Overview

### For Single Run
```
[0] (Recommended) Download Dorado model(s) once → persistent model directory
    ↓
POD5 files (single run)
    ↓
[2] run_dorado_basecaller.sh → Basecalled BAM (all barcodes)
    ↓
[4] run_dorado_demux.sh → Per-barcode FASTQ.gz files
```

*(If basecalling times out or fails, use Step 3: `resume_dorado_basecaller.sh`.)*

### For Multiple Runs (Recommended)
```
[0] (Recommended) Download Dorado model(s) once → persistent model directory
    ↓
POD5 files (multiple runs mixed)
    ↓
[1] Organize by run ID → Separate directories per run
    ↓
[2] run_dorado_basecaller.sh → One BAM per run
    ↓
[4] run_dorado_demux.sh → Per-barcode FASTQ per run
    ↓
[5] Merge same barcodes across runs (optional)
```

---

## Configuration

Before running, edit each script to replace placeholders.

### Common Placeholders

- `<PROJECT>`: Your NCI project code (e.g., `xe2`)
- `<DORADO_INSTALLATION_PATH>`: Path to Dorado installation directory
  - Example: `/g/data/xe2/adria/apps/dorado-1.1.0-linux-x64`

**Model download / storage**
- `<DORADO_MODEL_DIR>`: Directory where you want to store Dorado models (persistent location recommended)
  - Example: `/g/data/xe2/adria/apps/dorado_persistent_models`
- `<DORADO_MODEL_NAME>`: Dorado model identifier (the exact string used by `dorado download --model`)
  - Example: `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`

**If your basecalling scripts use a model *path* (rather than a model name):**
- `<DORADO_MODEL_PATH>`: Path to the downloaded model directory
  - Example: `<DORADO_MODEL_DIR>/<DORADO_MODEL_NAME>`
  - (After running the download script, Dorado will create a subdirectory under `<DORADO_MODEL_DIR>` matching the model name.)

- `<BARCODE_KIT>`: Barcoding kit used for library preparation
  - Example: `SQK-NBD114-96` (96-sample barcoding kit)

### Script-Specific Placeholders

**dorado_model_download.sh**:
- `<DORADO_INSTALLATION_PATH>`
- `<DORADO_MODEL_DIR>`
- `<DORADO_MODEL_NAME>`

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

## Step 0: Download Dorado Model(s) (Recommended)

Download the model once into a persistent directory (e.g., `gdata`) so basecalling jobs don’t re-download models and can start faster.

```bash
# Edit script to set paths
nano dorado_model_download.sh

# Run download (no GPU required)
bash dorado_model_download.sh

# Check models were downloaded
ls -lh <DORADO_MODEL_DIR>
```

**Notes**:
- If downloads fail due to network restrictions on your system, download the model on a machine with internet access and copy the model directory into `<DORADO_MODEL_DIR>`.
- If your basecalling scripts support it, you can also export the models directory so Dorado can find models automatically:
  ```bash
  export DORADO_MODELS_DIRECTORY="<DORADO_MODEL_DIR>"
  ```

---

## Important: Organizing Multiple Sequencing Runs

### Best Practice: Process Each Run Separately

If you have POD5 files from multiple sequencing runs in the same directory, **organize them by run BEFORE basecalling**.

#### Identifying Different Runs

POD5 files contain run IDs in their filenames. For example:
```
PBE46605_skip_0df7786f_bf504bcd_0.pod5    # Run 1
PBE46605_skip_ba040519_dfe5b02b_0.pod5    # Run 2
PBE46605_skip_c8f27651_3dbe505c_0.pod5    # Run 3
```

The format is: `<flowcell>_<mode>_<run_id>_<unique_id>_<number>.pod5`

## Step 1: Organize POD5 Files by Run
```bash
#!/bin/bash
# Organize POD5 files into run-specific directories

SOURCE_DIR="/path/to/mixed/pod5/files"
OUTPUT_BASE="/path/to/organized/pod5"

mkdir -p "$OUTPUT_BASE"

# Extract unique run IDs
cd "$SOURCE_DIR"
for pod5 in *.pod5; do
    # Extract run ID (third underscore-separated field)
    run_id=$(echo "$pod5" | cut -d'_' -f3)

    # Create directory for this run
    run_dir="$OUTPUT_BASE/run_${run_id}"
    mkdir -p "$run_dir"

    # Move or copy file
    mv "$pod5" "$run_dir/"  # Use 'cp' if you want to keep originals
done

echo "POD5 files organized by run:"
ls -d "$OUTPUT_BASE"/run_*
```

#### Expected Structure After Organization
```
organized/
├── run_0df7786f/
│   ├── PBE46605_skip_0df7786f_bf504bcd_0.pod5
│   ├── PBE46605_skip_0df7786f_bf504bcd_1.pod5
│   └── ...
├── run_ba040519/
│   ├── PBE46605_skip_ba040519_dfe5b02b_0.pod5
│   └── ...
└── run_c8f27651/
    ├── PBE46605_skip_c8f27651_3dbe505c_0.pod5
    └── ...
```

---

## Step 2: Initial Basecalling

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

## Step 3: Resume Incomplete Basecalling (If Needed)

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

## Step 4: Demultiplexing and FASTQ Conversion

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

### Common (Simplex) Models (R10.4.1 chemistry)

| Model                                      | Accuracy | Speed     | Use Case              |
|--------------------------------------------|----------|-----------|-----------------------|
| dna_r10.4.1_e8.2_400bps_sup@v5.2.0        | Highest  | Slow      | Publication-quality   |
| dna_r10.4.1_e8.2_400bps_hac@v5.2.0        | High     | Medium    | Standard analysis     |
| dna_r10.4.1_e8.2_400bps_fast@v5.2.0       | Moderate | Fast      | Quick QC              |

**sup** (super-accurate) is recommended for metagenomics to minimize assembly errors.

### Modified-Basecalling Models (Examples)

These are model identifiers you can use directly with `dorado_model_download.sh`:

- `dna_r10.4.1_e8.2_400bps_hac@v5.2.0_6mA@v1`
- `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2`

---

## Troubleshooting

### Job Killed or Timed Out

**Symptoms**: Job stops before completion, incomplete BAM file

**Solution**: Use `resume_dorado_basecaller.sh`
1. Note the incomplete BAM file path
2. Edit resume script with paths
3. Submit resume job
4. Repeat if necessary

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
2. Pre-download the model into a persistent directory using `dorado_model_download.sh`
3. Ensure Dorado can find the downloaded model (e.g., `export DORADO_MODELS_DIRECTORY=<DORADO_MODEL_DIR>`)
4. Check POD5 file integrity
5. Increase `--chunk-size` if enough GPU memory

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
# 0. (Recommended) Download Dorado model(s) once
bash scripts/00-basecalling/dorado_model_download.sh

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
# Step 0: Download model(s) once
bash dorado_model_download.sh

# Step 2: Initial basecalling (48 hours on V100)
# (Skip Step 1 if your POD5 are already in a single-run directory)
qsub run_dorado_basecaller.sh

# Check job status
qstat -u $USER

# If job completes successfully:
ls <OUTPUT_BAM_DIR>/*.bam

# Step 4: Demultiplex (4-8 hours)
qsub run_dorado_demux.sh

# Verify outputs
ls <OUTPUT_FASTQ_DIR>/*.fastq.gz | wc -l  # Should match number of barcodes used

# Proceed to quality control
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
