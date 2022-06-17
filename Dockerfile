FROM rocker/verse:4.2.0

# Update apt-get and install other libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libbz2-dev \
    libgdal-dev \
    libgeos-dev \
    libglpk40 \
    liblzma-dev \
    libmagick++-dev \
    libproj-dev \
    libudunits2-dev \
    libxt-dev \
    python3-pip \
    python3-dev

# Install pyrefinebio v0.4.9
RUN pip3 install pyrefinebio==0.4.9

# R Bioconductor packages
RUN Rscript -e "options(warn = 2); BiocManager::install(c( \
    'Biobase', \
    'BiocStyle', \
    'EnsDb.Hsapiens.v86', \
    'ensembldb', \
    'leukemiasEset', \
    'switchBox'), \
    update = FALSE, \
    version = 3.15)"

# R packages
RUN install2.r --error --deps TRUE \
    caret \
    here \
    multiclassPairs \
    optparse

