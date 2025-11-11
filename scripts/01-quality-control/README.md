# Quality Control Scripts

PBS scripts for Oxford Nanopore read quality assessment and filtering.

## Prerequisites

- Conda environment with: `nanoplot`, `fastqc`, `multiqc`, `chopper`
- PBS/Torque cluster with jobfs support

## Configuration

Before running, edit each script to replace placeholders:

- `<PROJECT>`: Your NCI project code (e.g., `xe2`)
- `<CONDA_PATH>`: Path to your conda installation
- `<QC_ENV_PATH>`: Path to your QC conda environment
- `<INPUT_FASTQ_DIR>`: Directory containing input FASTQ files
- `<OUTPUT_BASE_DIR>`: Base directory for results
- `<OUTPUT_FILTERED_DIR>`: Directory for filtered reads (chopper only)

## Workflow

### Step 1: NanoPlot (Read Statistics)
```bash
# Edit script to set paths
nano run_nanoplot.sh

# Submit job
qsub run_nanoplot.sh
```

**Output**: `<OUTPUT_BASE_DIR>/nanoplot/<sample>.txt` and `<sample>.html`

---

### Step 2: FastQC (Per-base Quality)
```bash
# Edit script to set paths
nano run_fastqc.sh

# Submit job
qsub run_fastqc.sh
```

**Output**: `<OUTPUT_BASE_DIR>/fastqc/<sample>_fastqc.html` and `.zip`

---

### Step 3: Chopper (Quality Filtering)
```bash
# Edit script to set paths and filtering parameters
nano filter_chopper.sh

# Adjust parameters:
MIN_LENGTH=2000    # Minimum read length (bp)
MIN_QUALITY=18     # Minimum average quality score

# Submit job
qsub filter_chopper.sh
```

**Output**: `<OUTPUT_FILTERED_DIR>/<sample>_filtered.fastq.gz`

---

### Step 4: MultiQC (Aggregate Report)

After running NanoPlot and/or FastQC:
```bash
# Edit script to set paths
nano multiqc_report.sh

# Submit job
qsub multiqc_report.sh
```

**Output**: `<OUTPUT_BASE_DIR>/multiqc/multiqc_report.html`

## Example: Complete QC Workflow
```bash
# 1. Run NanoPlot on raw reads
qsub run_nanoplot.sh

# 2. Run FastQC on raw reads
qsub run_fastqc.sh

# 3. Wait for jobs to complete, then aggregate
qsub multiqc_report.sh

# 4. Filter reads based on QC results
qsub filter_chopper.sh

# 5. Run QC on filtered reads
# (Edit scripts to point to filtered directory, then rerun)
qsub run_nanoplot.sh
qsub run_fastqc.sh
qsub multiqc_report.sh
```

## Resource Requirements

| Script          | Memory | CPUs | Walltime | Jobfs |
|-----------------|--------|------|----------|-------|
| run_nanoplot    | 64 GB  | 8    | 5 hours  | 400 GB|
| run_fastqc      | 64 GB  | 8    | 5 hours  | 400 GB|
| filter_chopper  | 32 GB  | 8    | 10 hours | 400 GB|
| multiqc_report  | 16 GB  | 1    | 1 hour   | N/A   |

## Troubleshooting

### FastQC Java Heap Space Error

If FastQC crashes with `OutOfMemoryError`:
- Increase `_JAVA_OPTIONS` in `run_fastqc.sh`
- Currently set to `-Xmx48G` (48 GB)

### NanoPlot Empty File Error

If NanoPlot fails with "no reads found":
- Check input file format (gzipped FASTQ expected)
- Verify reads start with `@` character
- Try uncompressed FASTQ

### MultiQC No Data

If MultiQC report is empty:
- Verify FastQC/NanoPlot directories exist
- Check file permissions
- Ensure FastQC/NanoPlot jobs completed successfully

## Output Structure
```
<OUTPUT_BASE_DIR>/
├── nanoplot/
│   ├── sample1.txt
│   ├── sample1.html
│   ├── sample2.txt
│   └── sample2.html
├── fastqc/
│   ├── sample1_fastqc.html
│   ├── sample1_fastqc.zip
│   ├── sample2_fastqc.html
│   └── sample2_fastqc.zip
└── multiqc/
    ├── multiqc_report.html
    └── multiqc_data/

<OUTPUT_FILTERED_DIR>/
├── sample1_filtered.fastq.gz
└── sample2_filtered.fastq.gz
```

## Notes

- Scripts use jobfs for temporary storage (faster I/O)
- Input files are copied to jobfs before processing
- Results are copied back to final output directory
- Temporary files are cleaned up automatically
- All scripts include error handling (`set -euo pipefail`)