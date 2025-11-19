# Long-Read Soil Metagenomics Pipeline

Workflow for high-quality genome assembly and differential abundance analysis 
from Oxford Nanopore soil metagenomic data.

## Overview

This repository contains the complete bioinformatic pipeline used in:

> Wu, H. (2025). Long-Read Soil Metagenomic Assembly and Differential 
> Abundance Analysis Across Plant Shoot Biomass Gradients. *Master's Thesis*, 
> Australian National University.

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