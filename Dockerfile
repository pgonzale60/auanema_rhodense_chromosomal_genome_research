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
    libjpeg-dev \
    gfortran \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install corporate SSL certificate if present (for corporate network)
COPY zscaler.crt /usr/local/share/ca-certificates/zscaler.crt
RUN update-ca-certificates && \
    echo "=== Certificate installed ===" && \
    ls -lh /etc/ssl/certs/ca-certificates.crt

# Set SSL environment variables for R to use system certificates
ENV CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
ENV R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE

# Configure R to use system certificates via Renviron.site
# These settings ensure renv uses the correct certificate bundle
RUN echo "CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt" > /usr/local/lib/R/etc/Renviron.site && \
    echo "SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt" >> /usr/local/lib/R/etc/Renviron.site && \
    echo "R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE" >> /usr/local/lib/R/etc/Renviron.site && \
    echo "=== R environment configuration ===" && \
    cat /usr/local/lib/R/etc/Renviron.site

# Verify R can see the SSL settings before attempting package installation
RUN echo "=== Testing R SSL configuration ===" && \
    R --vanilla --quiet -e "cat('CURL_CA_BUNDLE:', Sys.getenv('CURL_CA_BUNDLE'), '\n')" && \
    R --vanilla --quiet -e "cat('SSL_CERT_FILE:', Sys.getenv('SSL_CERT_FILE'), '\n')" && \
    R --vanilla --quiet -e "cat('R_LIBCURL_SSL_REVOKE_BEST_EFFORT:', Sys.getenv('R_LIBCURL_SSL_REVOKE_BEST_EFFORT'), '\n')"

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
