# Long-Read Soil Metagenomics Pipeline

Workflow for high-quality genome assembly and differential abundance analysis 
from Oxford Nanopore soil metagenomic data.

## Overview

This repository contains the complete bioinformatic pipeline used in:

> Wu, H. (2025). Long-Read Soil Metagenomic Assembly and Differential 
> Abundance Analysis Across Plant Shoot Biomass Gradients. *Master's Thesis*, 
> Australian National University.

**📄 [Read the full thesis (PDF)](https://github.com/AdriaWu/soil-metagenome-workflow/blob/main/docs/HeyueWu_masters_thesis.pdf)**

## Data Availability

All raw sequencing data from this study are publicly available through NCBI:

**🔗 [NCBI BioProject PRJNA1400213](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1400213)**

- **358 BioSamples** across three experimental datasets
- **193 Gbp** of Oxford Nanopore R10 sequencing data
- Quality-filtered reads (≥300 bp, mean Q≥15)
- Comprehensive sample metadata included

### Datasets

The BioProject encompasses three complementary experiments:

1. **Agricultural soil time series** - Boorowa (NSW) soil cores sampled in 2019 and 2023 with soil chemistry measurements
2. **Plant productivity gradient** - Rhizosphere microbiomes across shoot biomass gradients in multiple crop hosts
3. **Forest vs. river soil comparison** - ACT soils under different fertilizer treatments

### Key Features

- **High-contiguity assembly** using nanoMDBG optimized for ONT R10 data
- **Rigorous contamination removal** through competitive mapping
- **Compositionally-aware differential abundance testing** with empirical FDR
- **Genome-resolved ecological interpretation** linking taxa to functional traits

### Workflow Summary
```
Raw ONT Reads
    ↓
Quality Control & Filtering (Chopper, NanoPlot)
    ↓
Contamination Removal (minimap2 + competitive classification)
    ↓
Co-assembly (nanoMDBG)
    ↓
Assembly QC (MetaQUAST)
    ↓
Read Mapping & Count Matrix
    ↓
Differential Abundance (edgeR + empirical FDR)
    ↓
Taxonomic Annotation & Visualization
```

## Quick Start

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- 32+ GB RAM recommended
- ~500 GB disk space for intermediate files

### Installation
```bash
# Clone repository
git clone https://github.com/AdriaWu/soil-metagenome-workflow.git
cd soil-metagenome-workflow

# Create conda environment
conda env create -f environment.yml
conda activate soil-meta
```

## Key Results

- **556,967 contigs** totaling 7.69 Gbp
- **86 contigs >1 Mb** (longest: 4.52 Mb)
- **923 circular replicons** (complete genetic elements)
- **20% higher N50** vs. global SMAG catalogue
- **194 differentially abundant contigs** across biomass gradient

## Citation

**📄 [Download full thesis (PDF)](https://github.com/AdriaWu/soil-metagenome-workflow/blob/main/docs/HeyueWu_masters_thesis.pdf)**

If you use this workflow, please cite:
```bibtex
@mastersthesis{wu2025longread,
  author = {Wu, Heyue},
  title = {Long-Read Soil Metagenomic Assembly and Differential Abundance 
           Analysis Across Plant Shoot Biomass Gradients},
  school = {Australian National University},
  year = {2025},
  type = {Master's Thesis}
}
```

## License

[MIT License](LICENSE)

## Contact

- **Author**: Heyue Wu
- **Supervisor**: Professor Justin Borevitz
- **Institution**: Research School of Biology, Australian National University

## Acknowledgments

Special thanks to the Borevitz Lab, particularly Viraj Kolhapuri, Alek Meade, 
Saige Waugh, Ashley Jones, Christopher Bradley for their contributions.