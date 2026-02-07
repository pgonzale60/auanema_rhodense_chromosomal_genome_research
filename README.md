# Chromosomal Genome Analysis Repository

Reproducible analysis pipeline for genome characterization, including chromatin features, tandem repeats, and gene annotation.

## Overview

This repository contains the complete analysis workflow for a chromosomal genome project, with a focus on:
- **Genome characterization** (assembly QC, structural features)
- **Chromatin organization** (Hi-C analysis, 3D structure)
- **Tandem repeats** (TRF family identification and classification)
- **Gene annotation** (functional annotation, orthology)

## Core Analysis Scripts

The following scripts generate the main paper figures:

| Script | Output | Purpose |
|--------|--------|---------|
| `scripts/fig_2_hiC_and_genomic_features.R` | `report/figures/fig_2_hiC_and_genomic_features.pdf` | Figure 2: Hi-C + integrated genomic features |
| `scripts/fig_3_tandem_repeats.R` | `report/figures/fig_3_tandem_repeats.pdf` | Figure 3: Tandem repeat family distribution |
| `scripts/fig_4_motif_analysis.R` | `report/figures/fig_4_motif_analysis.pdf` | Figure 4: Motif analysis visualization |
| `scripts/summarize_gene_annotation.R` | Console output / summary tables | Gene annotation statistics |

## Data Dependencies

### Genome Assembly
The primary genome assembly is **not included** in this repository.
- **NCBI BioProject:** PRJEB75808
- Download `nxAuaRhod1_1.primary.fa.gz` and place in repository root (gitignored)

### Hi-C Data
- **File:** `analyses/hi-C/cooltools/nxAuaRhod1_1.mcool` (53MB, tracked in repository)
- Multi-resolution contact matrix for chromatin organization analysis

### Annotation Data
All compressed annotation files are tracked in `analyses/`:
- Gene models: `analyses/genes/GTFs/nxAuaRhod1_1.gff.gz`
- Orthology data: `analyses/genes/orthofinder/`
- Functional annotations: `analyses/genes/annotation/`
- Repeat analysis: `analyses/genome_features/repeats/TRF/`
- Motif analysis: `analyses/diminution/`

## Quick Start

### GitHub Codespaces / VS Code Dev Container (Recommended)

The repository ships with a `.devcontainer` that provisions:
- R 4.4.2 (rocker/r-ver)
- System libraries for Bioconductor and tidyverse
- `renv` + cache support (packages installed on demand)
- VS Code R tooling

**Codespaces**
1. Click the green "Code" button on GitHub → Codespaces tab
2. Create a new Codespace on `main`
3. Wait ~3–4 minutes for the container image to build

**Local VS Code**
1. Install Docker Desktop and the "Dev Containers" extension
2. Open the repo in VS Code
3. Run "Dev Containers: Reopen in Container" from the command palette

Once inside the container, open R and restore packages:
```r
renv::restore()
```

The first restore may take a few minutes as Bioconductor packages are compiled; subsequent runs are faster.

## Local R Setup (Without Docker)

If you prefer to run locally:

1. Install R 4.4.2
2. Install system dependencies (example for macOS):
   ```bash
   brew install gfortran
   ```
3. Open R in this directory:
   ```r
   renv::restore()  # Install all packages from renv.lock
   ```

**Note:** You may encounter the gfortran library path issue on macOS. Docker is the recommended approach.

## Reproducibility

This project uses [renv](https://rstudio.github.io/renv/) to ensure reproducible R package versions:
- R version: 4.4.2
- Bioconductor: 3.20
- All package versions locked in `renv.lock`

The Docker setup ensures system-level reproducibility (gfortran, system libraries, etc.).

## Citation

If you use this code or data in your research, please cite this repository:

```
Gonzalez de la Rosa, P. (2025). Chromosomal Genome Analysis.
GitHub repository: https://github.com/pgonzale60/auanema_rhodense_chromosomal_genome_research
```

A peer-reviewed publication is in preparation.

## Contact

For questions or collaboration inquiries:
- **Pablo Gonzalez de la Rosa**: pgonzale60@gmail.com | [@pgonzale60](https://github.com/pgonzale60)
- **Mark Blaxter**: mb35@sanger.ac.uk

## License

This project is dual-licensed following Wellcome Sanger Institute practices:
- **Code** (R/Python scripts): AGPL-3.0
- **Data** (genome annotations, results): CC BY 4.0

See the [LICENSE](LICENSE) file for full details.

Copyright (c) 2025 Genome Research Ltd. (Wellcome Sanger Institute)
