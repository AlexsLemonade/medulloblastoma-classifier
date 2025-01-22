FROM bioconductor/bioconductor_docker:3.16

# set a name for the conda environment
ARG ENV_NAME=medulloblastoma-classifier

# set environment variables to install conda
ENV PATH="/opt/conda/bin:${PATH}"

# Install conda via miniforge
# adapted from https://github.com/conda-forge/miniforge-images/blob/master/ubuntu/Dockerfile
RUN curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" -o /tmp/miniforge.sh \
  && bash /tmp/miniforge.sh -b -p /opt/conda \
  && rm -f /tmp/miniforge.sh \
  && conda clean --tarballs --index-cache --packages --yes \
  && find /opt/conda -follow -type f -name '*.a' -delete \
  && find /opt/conda -follow -type f -name '*.pyc' -delete \
  && conda clean --force-pkgs-dirs --all --yes

# Activate conda environments in bash
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
  && echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/skel/.bashrc \
  && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# Install conda-lock
RUN conda install --channel=conda-forge --name=base conda-lock \
  && conda clean --all --yes

# Copy conda lock file to image
COPY conda-lock.yml conda-lock.yml

# restore from conda-lock.yml file and clean up to reduce image size
RUN conda-lock install -n ${ENV_NAME} conda-lock.yml \
  && conda clean --all --yes

# Activate conda environment on bash launch
RUN echo "conda activate ${ENV_NAME}" >> ~/.bashrc

# Use renv for R packages
WORKDIR /usr/local/renv
ENV RENV_CONFIG_CACHE_ENABLED=FALSE
RUN Rscript -e "install.packages('renv')"

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

# Copy the renv.lock file from the host environment to the image
COPY renv.lock renv.lock

# restore from renv.lock file and clean up to reduce image size
RUN Rscript -e 'renv::restore()' \
  && rm -rf ~/.cache/R/renv \
  && rm -rf /tmp/downloaded_packages \
  && rm -rf /tmp/Rtmp*

WORKDIR /home/rstudio
