# Auanema rhodense Genome Analysis

Reproducible analysis of *Auanema rhodense* gene annotation, orthology, and genomic features.

## Overview

This repository contains R scripts and data for analyzing the *Auanema rhodense* genome (assembly nxAuaRhod1.1), with a focus on:
- Gene annotation (InterProScan, EggNOG-mapper)
- Orthology analysis (OrthoFinder)
- Gene expression (STAR/StringTie)
- Transposon-related genes
- Genome-wide features

## Quick Start with Docker

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

If you use this code or data, please cite:

[Your publication details here]

## Contact

[Your contact information]

## License

[Add license here]
