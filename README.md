# Auanema rhodense Genome Analysis

Reproducible analysis of *Auanema rhodense* gene annotation, orthology, and genomic features.

## Overview

This repository contains R scripts and data for analyzing the *Auanema rhodense* genome (assembly nxAuaRhod1.1), with a focus on:
- Gene annotation (InterProScan, EggNOG-mapper)
- Orthology analysis (OrthoFinder)
- Gene expression (STAR/StringTie)
- Transposon-related genes
- Genome-wide features

## Quick Start Options

### Option 1: GitHub Codespaces (Recommended)

The easiest way to get started is using GitHub Codespaces, which provides a cloud-based development environment:

1. Click the green "Code" button on the GitHub repository
2. Select "Codespaces" tab
3. Click "Create codespace on main"

The environment will automatically set up with:
- R 4.4.2
- RStudio Server (accessible via browser)
- All system dependencies (gfortran, GSL, etc.)
- renv for package management
- VSCode R extensions

After the container starts (2-3 minutes), install R packages as needed:
```r
# Install all packages from renv.lock (takes ~15-20 minutes first time)
renv::restore()

# Or install specific packages only
renv::install("tidyverse")
renv::install("DESeq2")
```

### Option 2: Local VSCode Dev Container

If you have Docker and VSCode installed locally:

1. Install the "Dev Containers" extension in VSCode
2. Open the repository folder in VSCode
3. Press `Cmd+Shift+P` (Mac) or `Ctrl+Shift+P` (Windows/Linux)
4. Select "Dev Containers: Reopen in Container"

Container builds in 2-3 minutes. Then run `renv::restore()` in R to install packages.

### Option 3: Docker Compose (Local)

### Prerequisites
- [Docker](https://www.docker.com/get-started) installed on your system
- [Docker Compose](https://docs.docker.com/compose/install/) (usually included with Docker Desktop)
- 4GB+ RAM available for Docker

### 1. Clone the Repository
```bash
git clone https://github.com/pgonzale60/auanema_rhodense_chromosomal_genome_research.git
cd auanema_rhodense_chromosomal_genome_research
```

### 2. Build the Docker Image
```bash
docker-compose build
```

This will:
- Create an R 4.4.2 environment
- Install all system dependencies (including gfortran for Matrix package)
- Restore all R packages from `renv.lock` (Bioconductor 3.20 + tidyverse)
- Takes ~10-15 minutes on first build

### 3. Run R in the Container
```bash
docker-compose run --rm r-analysis R
```

You're now in an R session with all packages installed!

### 4. Run an Analysis Script
```bash
docker-compose run --rm r-analysis Rscript scripts/summarize_gene_annotation.R
```

## Alternative: Interactive R Session

Start a bash shell in the container:
```bash
docker-compose run --rm r-analysis bash
```

Then run R or any scripts:
```bash
R
# or
Rscript scripts/summarize_gene_annotation.R
```

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
├── Dockerfile                        # Docker build instructions
└── docker-compose.yml                # Docker orchestration
```

## Data Files

All data files are compressed (gzip) to reduce repository size. R's `read_tsv()` automatically handles `.gz` files.

Key data files:
- `analyses/genes/annotation/nxAuaRhod1_1.iprscan.guoying.tsv.gz` - InterProScan results
- `analyses/genes/annotation/nxAuaRhod1_1.eggnog_mapper.tsv.gz` - EggNOG-mapper annotations
- `analyses/genes/orthofinder/Orthogroups/Orthogroups.tsv.gz` - Orthogroup assignments
- `analyses/genes/star_stringtie/arhod_star_stringtie_gene_expression.tsv.gz` - Gene expression

## Known Issues

### Missing `metadata` Object

**Issue:** The `summarize_gene_annotation.R` script references a `metadata` object on line 29 that is not defined in the script.

```r
left_join(select(metadata, sample, name), by = c("SRA" = "sample"))
```

**Workaround:** You'll need to either:
1. Load the metadata from an external file (if available)
2. Comment out line 29 if gene expression metadata is not needed for your analysis
3. Create a minimal metadata object with the required columns

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
