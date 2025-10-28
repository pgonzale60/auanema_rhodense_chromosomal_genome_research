# Auanema rhodense Genome Analysis

Reproducible analysis of *Auanema rhodense* gene annotation, orthology, and genomic features.

## Overview

This repository contains R scripts and data for analyzing the *Auanema rhodense* genome (assembly nxAuaRhod1.1), with a focus on:
- Gene annotation (InterProScan, EggNOG-mapper)
- Orthology analysis (OrthoFinder)
- Gene expression (STAR/StringTie)
- Transposon-related genes
- Genome-wide features

## Quick Start

### GitHub Codespaces / VS Code Dev Container (Recommended)

The repository ships with a `.devcontainer` that provisions:
- R 4.4.2 (rocker/r-ver)
- System libraries required by Bioconductor and tidyverse
- `renv` + cache support (packages installed on demand)
- VS Code R tooling (languageserver, debugger, session watcher)

**Codespaces**
1. Click the green “Code” button on GitHub → Codespaces tab.
2. Create a new Codespace on `main`.
3. Wait ~3–4 minutes for the container image to build.

**Local VS Code**
1. Install Docker Desktop and the “Dev Containers” extension.
2. Open the repo in VS Code.
3. Run “Dev Containers: Reopen in Container” from the command palette.

Once inside the container (cloud or local), open R and hydrate packages when needed:
```r
renv::restore()
```

The first restore compiles numerous Bioconductor packages; subsequent runs are faster thanks to the shared cache.

## Repository Structure

```
.
├── scripts/                    # R analysis scripts
│   ├── summarize_gene_annotation.R    # Main gene annotation analysis
│   ├── orthologs.R                    # Orthology analysis
│   └── ...
├── analyses/
│   ├── genes/
│   │   ├── annotation/               # InterProScan, EggNOG results
│   │   ├── orthofinder/             # OrthoFinder orthogroups
│   │   └── star_stringtie/          # Gene expression data
├── raw/                              # Reference data files
├── renv.lock                         # R package versions (DO NOT MODIFY)
├── .Rprofile                         # Activates renv
├── .devcontainer/                    # Dev container definition (Dockerfile + settings)
```

## Data Files

All data files are compressed (gzip) to reduce repository size. R's `read_tsv()` automatically handles `.gz` files.

Key data files:
- `analyses/genes/annotation/nxAuaRhod1_1.iprscan.guoying.tsv.gz` - InterProScan results
- `analyses/genes/annotation/nxAuaRhod1_1.eggnog_mapper.tsv.gz` - EggNOG-mapper annotations
- `analyses/genes/orthofinder/Orthogroups/Orthogroups.tsv.gz` - Orthogroup assignments
- `analyses/genes/star_stringtie/arhod_star_stringtie_gene_expression.tsv.gz` - Gene expression

## Reproducibility

This project uses [renv](https://rstudio.github.io/renv/) to ensure reproducible R package versions:
- R version: 4.4.2
- Bioconductor: 3.20
- All package versions locked in `renv.lock`

The Docker setup ensures system-level reproducibility (gfortran, system libraries, etc.).

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

## Citation

If you use this code or data in your research, please cite this repository:

```
Gonzalez de la Rosa, P. (2025). Auanema rhodense Genome Analysis.
GitHub repository: https://github.com/pgonzale60/auanema_rhodense_chromosomal_genome_research
```

A peer-reviewed publication is in preparation.

## Contact

For questions or collaboration inquiries:
- **Pablo Gonzalez de la Rosa**: pgonzale60@gmail.com | [@pgonzale60](https://github.com/pgonzale60)
- **Mark Blaxter**: mb35@sanger.ac.uk

## License

This project is dual-licensed following Wellcome Sanger Institute practices:
- **Code** (R scripts): AGPL-3.0
- **Data** (genome annotations, results): CC BY 4.0

See the [LICENSE](LICENSE) file for full details.

Copyright (c) 2025 Genome Research Ltd. (Wellcome Sanger Institute)
