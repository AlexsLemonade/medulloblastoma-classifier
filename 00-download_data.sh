# Download data (metadata and expression data) associated with several projects
#
# Steven Foltz
# November 2022
#
# Usage: 00_download_data.sh --data_sources_file data/data_sources.tsv
# where data/data_sources.tsv is a TSV file with accession, data_source, and url columns

#!/bin/bash
set -euo pipefail

# change to the directory of this script
cd "$(dirname "${BASH_SOURCE[0]}")"

# set data directory
data="data"

# set input data sources file
data_sources_file="data/data_sources.tsv"

# allow for alternative input data sources file
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

if [[ ! -f $data_sources_file ]]; then

  echo Input data_sources_file $data_sources_file does not exist.
  exit 1
  
fi

# read in each column of data sources file one line at a time
# data_source defines the action taken for each line
while read accession data_source url; do

  # skip over any header lines starting with hash
  [[ $accession =~ ^#.* ]] && continue

  if [[ $accession == "E-MTAB-292" ]]; then
  
    # E-MTAB-292 is part of ArrayExpress.
    # We can come back to this dataset in the future.
    # For now it is not clear how to easily download processed expression data.
    echo Skipping over E-MTAB-292.

  elif [[ $download_source == refine.bio ]]; then

    if [[ -d $data/$accession ]]; then
  
      echo Data for $accession already exists and will not be downloaded.
  
    else
    
      echo Downloading $accession from refine.bio...
      refinebio create-token -s
    
      refinebio download-dataset \
        --email-address steven.foltz@ccdatalab.org \
        --path ${data}/${accession}.zip \
        --experiments $accession \
        --aggregation EXPERIMENT \
        --transformation NONE \
        --skip-quantile-normalization True

      unzip -d ${data}/${accession} ${data}/${accession}.zip && rm -f ${data}/${accession}.zip
    
      sleep 1
      
    fi

  elif [[ $download_source == url ]]; then

    mkdir -p $data/$accession
    url_basename=$(basename $url)
    file_name=$data/$accession/$url_basename
    
    if [[ -f $file_name ]]; then
      
      echo File $file_name already exists and will not be downloaded.
  
    else
    
      echo Downloading $accession $url_basename...
      
      curl -o $file_name --silent $url
      
      sleep 1
      
    fi
  
  else
  
    echo Error with: $accession $download_source $url
    exit 1
    
  fi
  
done < $data_sources_file

################################################################################
# md5sum of downloaded data
################################################################################

# As needed, recreate data/md5_check_sums.tsv by running this command:
# find data -type f | grep -v "json\|md5" | xargs md5sum > data/md5_check_sums.tsv

# check md5 sums of downloaded files
echo Checking md5 sums of downloaded files...
md5sum --check --quiet data/md5_check_sums.tsv && echo All good!
