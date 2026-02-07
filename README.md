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
| `scripts/create_diminution_figure.R` | `report/figures/diminution_motif_figure.pdf` | Motif analysis visualization |
| `scripts/plot_trf_families_chromosomes.R` | `report/figures/trf_families_grs_winners.pdf` | Repeat family chromosomal distribution |
| `scripts/plot_hic_and_genomic_features.R` | `report/figures/fig_2_hiC_and_genomic_features.pdf` | Hi-C + integrated genomic features |
| `scripts/summarize_gene_annotation.R` | Console output / summary tables | Gene annotation statistics |

## Data Dependencies

### Genome Assembly
The primary genome assembly is **not included** in this repository (obtain from NCBI/ENA):
- NCBI BioProject: PRJNA###
- Download: `nxAuaRhod1_1.primary.fa.gz`
- Add to root directory (will be ignored by `.gitignore`)

### Hi-C Data
- Cool file (.mcool): Not tracked; regenerate from raw sequencing data or obtain from the authors
- Reference: See Hi-C processing scripts in `analyses/`

### Annotation Data
All compressed annotation files (<10MB) are tracked in `analyses/`:
- Gene models: `analyses/genes/GTFs/nxAuaRhod1_1.gff.gz`
- Orthology data: `analyses/genes/orthofinder/`
- Functional annotations: `analyses/genes/annotation/`

## Repository Structure

```
.
├── scripts/                                    # Analysis and figure scripts
│   ├── create_diminution_figure.R              # Motif visualization
│   ├── plot_trf_families_chromosomes.R         # Repeat family figures
│   ├── plot_hic_and_genomic_features.R         # Integrated Hi-C figures
│   ├── summarize_gene_annotation.R             # Annotation statistics
│   ├── summarize_eliminated_genes.R            # Eliminated region analysis
│   ├── summarize_repeats.R                     # Repeat summary
│   ├── telomere_addition_positions.R           # Telomere analysis
│   ├── orthologs.R                             # Orthology analysis
│   ├── annotate_trf_families.py                # TRF family annotation pipeline
│   └── filter_trf_by_rna.py                    # TRF filtering utility
│
├── analyses/
│   ├── genes/
│   │   ├── GTFs/                               # Gene models
│   │   ├── annotation/                         # Functional annotations
│   │   ├── orthofinder/                        # Orthology results
│   │   └── star_stringtie/                     # Expression data
│   │
│   ├── genome_features/
│   │   ├── repeats/
│   │   │   ├── TRF/                            # Tandem repeat family data
│   │   │   ├── structural_rna*.bed             # RNA region annotations
│   │   │   └── modDotPlot/                     # (gitignored - intermediate)
│   │   ├── nemaChromQC/                        # QC analyses (BUSCO, GC, RED, telomeres)
│   │   └── elim_coords/                        # Eliminated regions
│   │
│   ├── diminution/                             # Motif analysis
│   │   ├── fimo_out/                           # FIMO results
│   │   ├── meme/                               # MEME motif outputs
│   │   ├── grs_visualization/                  # Visualization coordinates
│   │   ├── nxAuaRhod1_1.GRS.bed                # Region coordinates
│   │   └── nxAuaRhod1_1.core.bed               # Core region
│   │
│   └── curated_GRS_coords/                     # Curated region coordinates
│
├── report/
│   └── figures/                                # Final figure outputs
│       ├── fig_2_hiC_and_genomic_features.pdf
│       ├── trf_families_grs_winners.pdf
│       ├── diminution_motif_figure.pdf
│       ├── figure_2_legend.md
│       └── figure_trf_families_grs_winners_legend.md
│
├── renv.lock                                   # R package versions (locked)
├── .Rprofile                                   # renv activation
├── .devcontainer/                              # Docker dev environment
└── README.md                                   # This file
```

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

## File Notes

**Data files are compressed** (gzip) to reduce repository size. R's `read_tsv()` automatically handles `.gz` files.

**Key untracked files** (document in `.gitignore` / data README how to obtain):
- Genome assembly: `nxAuaRhod1_1.primary.fa.gz*`
- Hi-C data: `.mcool` files

**Gitignored directories** (intermediate outputs):
- `analyses/genome_features/repeats/modDotPlot/` — intermediate dot plot data

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
