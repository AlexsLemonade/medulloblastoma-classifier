# Download data (metadata and expression data) associated with several projects
#
# Steven Foltz
# November 2022

#!/bin/bash
set -euo pipefail

# change to the directory of this script
cd "$(dirname "${BASH_SOURCE[0]}")"

# set data directory
data="data"
mkdir -p $data

################################################################################
# GSE124814
################################################################################

# GSE124814 is a series containing several datasets.
# We can use the harmonized metadata from the series but will download each
# dataset independently from refine.bio, except E-MTAB-292.

if [[ -d $data/GSE124814 ]]; then
  
  echo Metadata for series GSE124814 already exists and will not be downloaded.

else

  mkdir -p $data/GSE124814
  GSE124814_sample_descriptions_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124814/suppl/GSE124814_sample_descriptions.xlsx"
  GSE124814_sample_descriptions_basename=$(basename $GSE124814_sample_descriptions_url)
  
  echo Downloading metadata for series GSE124814...
  curl -o $data/GSE124814/$GSE124814_sample_descriptions_basename --silent $GSE124814_sample_descriptions_url

fi

# get medulloblastoma studies identified by GSE124814
while read accession; do

  if [[ -d $data/$accession ]]; then
  
    echo Data for $accession already exists and will not be downloaded.

  elif [[ $accession == GSE* ]]; then
  
    echo Downloading $accession...
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
    
  elif [[ $accession == "E-MTAB-292" ]]; then
  
    # E-MTAB-292 is part of ArrayExpress.
    # We can come back to this dataset in the future.
    # For now it is not clear how to easily download processed expression data.
    echo Skipping over E-MTAB-292.
    
  else
  
    echo $accession does not start with GSE or equal E-MTAB-292...
    
  fi
  
done < $data/GSE124814_experiment_accessions.tsv

################################################################################
# GSE164677 -- not in refine.bio
################################################################################

if [[ -d $data/GSE164677 ]]; then

  echo Data for GSE164677 already exists and will not be downloaded.
  
else

  echo Downloading GSE164677...
  mkdir $data/GSE164677
  
  GSE164677_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164677/suppl/GSE164677_Asian_MB_RNA-seq.txt.gz"
  GSE164677_basename=$(basename $GSE164677_url)
  
  curl -o $data/GSE164677/$GSE164677_basename --silent $GSE164677_url
  
fi

################################################################################
# OpenPBTA
################################################################################

if [[ -d $data/OpenPBTA ]]; then

  echo Data for OpenPBTA already exists and will not be downloaded.

else

  echo Downloading OpenPBTA...
  mkdir $data/OpenPBTA
  
  OpenPBTA_release="release-v22-20220505"
  
  OpenPBTA_metadata_url="https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/data/$OpenPBTA_release/pbta-histologies.tsv"
  OpenPBTA_genex_stranded_url="https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/data/$OpenPBTA_release/pbta-gene-expression-rsem-tpm.stranded.rds"
  OpenPBTA_genex_polya_url="https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/data/$OpenPBTA_release/pbta-gene-expression-rsem-tpm.polya.rds"
  OpenPBTA_metadata_basename=$(basename $OpenPBTA_metadata_url)
  OpenPBTA_genex_stranded_basename=$(basename $OpenPBTA_genex_stranded_url)
  OpenPBTA_genex_polya_basename=$(basename $OpenPBTA_genex_polya_url)
  
  curl -o $data/OpenPBTA/$OpenPBTA_metadata_basename --silent $OpenPBTA_metadata_url
  curl -o $data/OpenPBTA/$OpenPBTA_genex_stranded_basename --silent $OpenPBTA_genex_stranded_url
  curl -o $data/OpenPBTA/$OpenPBTA_genex_polya_basename --silent $OpenPBTA_genex_polya_url

fi

################################################################################
# md5sum of downloaded data
################################################################################

# As needed, recreate data/md5_check_sums.tsv by running this command:
# find data -type f | grep -v "json\|md5" | xargs md5sum > data/md5_check_sums.tsv

# check md5 sums of downloaded files
echo Checking md5 sums of downloaded files...
md5sum --check --quiet data/md5_check_sums.tsv && echo All good!
