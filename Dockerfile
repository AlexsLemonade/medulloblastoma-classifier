FROM rocker/verse:4.2.2

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
    'AnnotationHub', \
    'Biobase', \
    'BiocStyle', \
    'GSVA', \
    'leukemiasEset', \
    'switchBox'), \
    update = FALSE, \
    version = 3.16)"

# R packages (R-forge)
RUN install2.r --error --deps TRUE --repos http://r-forge.r-project.org \
    estimate
    
# R packages (CRAN)
RUN install2.r --error --deps TRUE --repos http://cran.r-project.org \
    caret \
    doParallel \
    glmnet \
    here \
    MM2S \
    multiclassPairs \
    optparse \
    patchwork

# Threading issue with preprocessCore::normalize.quantiles
# https://support.bioconductor.org/p/122925/#124701
# https://github.com/bmbolstad/preprocessCore/issues/1#issuecomment-326756305
# put this last with force = TRUE to ensure it is properly installed
RUN Rscript -e "options(warn = 2); BiocManager::install( \
    'preprocessCore', \
    configure.args = '--disable-threading', \
    force = TRUE, \
    update = FALSE, \
    version = 3.16)"

