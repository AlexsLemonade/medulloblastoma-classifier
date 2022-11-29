# Download the raw data and metadata associated with project GSE119926 found at:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119926
#
# Chante Bethell
# November 2022
#
# Usage: bash download_GSE119926.sh

#!/bin/bash
set -euo pipefail

# change to the directory of this script
cd "$(dirname "${BASH_SOURCE[0]}")"

# set data directory
data_dir="data/GSE119926"
mkdir -p ${data_dir}

# download GSE119926 raw expression data
wget -O ${data_dir}/GSE119926_RAW.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE119926&format=file'

# download GSE119926 metadata
wget -O ${data_dir}/GSE119926_metadata.txt.gz 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119926/matrix/GSE119926_series_matrix.txt.gz'
