# Auanema rhodense genome analysis - Reproducible R environment
# Based on Rocker project: https://rocker-project.org/

FROM rocker/r-ver:4.4.2

# Install system dependencies for R packages
# These are needed for: Matrix, tidyverse, Bioconductor packages, XML parsing, etc.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg62-turbo-dev \
    gfortran \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Copy renv files first (for Docker layer caching)
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/.gitignore renv/.gitignore
COPY renv/settings.json renv/settings.json

# Install renv and restore packages
# This will install all packages from renv.lock
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))" \
    && R -e "renv::restore()"

# Copy the rest of the project
COPY . .

# Set default command to R
CMD ["R"]
