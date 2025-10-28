# Minimal dev container for Auanema rhodense analysis
# Uses rocker/r-ver (no RStudio) and leaves renv::restore() to run manually.

FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies commonly needed by Bioconductor / tidyverse packages
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
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libgsl-dev \
    gfortran \
    git \
    make \
    pkg-config \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Install renv and the language server for editor integration.
RUN R -e "install.packages(c('renv', 'languageserver'), repos = c(CRAN = 'https://cloud.r-project.org'))"

# Create a non-root user that matches Codespaces expectations
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} -m ${USERNAME} \
    && usermod -aG sudo ${USERNAME} \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Provide a consistent workdir for VS Code / Codespaces volumes
WORKDIR /workspace
RUN chown ${USERNAME}:${USERNAME} /workspace
