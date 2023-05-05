FROM rocker/r-ver:4.2.2

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


# Renv install all packages
WORKDIR /renv
COPY renv.lock renv.lock

RUN Rscript -e "install.packages(c('renv', 'BiocManager'))"
RUN Rscript -e "renv::restore()" && \
      rm -rf ~/.local/share/renv && \
      rm -rf /tmp/downloaded_packages && \
      rm -rf /tmp/Rtmp*

ENV RENV_DISABLED TRUE
WORKDIR /home
