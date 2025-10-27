# Auanema rhodense genome analysis - Dev Container
# Based on Rocker project: https://rocker-project.org/

FROM rocker/rstudio:4.4.2

# Install system dependencies for R packages
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
    libjpeg-dev \
    gfortran \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Give rstudio user sudo access for installing additional packages
RUN echo "rstudio ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Install renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Install R language server for VSCode integration (before renv::restore())
RUN R -e "install.packages('languageserver', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Set working directory
WORKDIR /workspace

# Copy renv configuration files
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/.gitignore renv/.gitignore
COPY renv/settings.json renv/settings.json

# Restore all R packages from renv.lock during build
# This caches packages in the Docker image layer
RUN R -e "renv::restore()"

# Change ownership to rstudio user
RUN chown -R rstudio:rstudio /workspace

# Default command (RStudio server for web interface)
CMD ["/init"]
