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
# dataset independently from refine.bio (and E-MTAB-292 from ArrayExpress).

if [[ -d $data/GSE124814 ]]; then
  
  echo Metadata for series GSE124814 already exists and will not be downloaded.

else

  mkdir -p $data/GSE124814
  GSE124814_sample_descriptions_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124814/suppl/GSE124814_sample_descriptions.xlsx"
  GSE124814_sample_descriptions_basename=$(basename $GSE124814_sample_descriptions_url)
  curl -o $data/GSE124814/$GSE124814_sample_descriptions_basename --silent $GSE124814_sample_descriptions_url
  
fi

# re-run metadata processing on GSE124814
Rscript -e "rmarkdown::render('GSE124814_metadata.Rmd')"

# get medulloblastoma studies identified by GSE124814
while read accession; do

  if [[ -d $data/$accession ]]; then
  
    echo Data for $accession already exists and will not be downloaded.

  elif [[ $accession == GSE* ]]; then
  
    echo Downloading $accession
    refinebio create-token -s
    refinebio download-dataset \
      --email-address steven.foltz@ccdatalab.org \
      --path ${data}/${accession}.zip \
      --experiments $accession \
      --aggregation EXPERIMENT \
      --transformation NONE \
      --skip-quantile-normalization True

    unzip -d ${data}/${accession} ${data}/${accession}.zip && rm -f ${data}/${accession}.zip
    
    sleep 5
    
  elif [[ $accession == "E-MTAB-292" ]]; then
  
    echo Downloading $accession
    mkdir -p $data/$accession
    idf_url="https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-292/E-MTAB-292.idf.txt"
    sdrf_url="https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-292/E-MTAB-292.sdrf.txt"
    genex_url="https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-292/E-MTAB-292.processed.1.zip/rma-gene-core.summary.txt"
    curl -o $data/$accession/$(basename $idf_url) --silent $idf_url
    curl -o $data/$accession/$(basename $sdrf_url) --silent $sdrf_url
    curl -o $data/$accession/$(basename $genex_url) --silent $genex_url
  
  else
  
    echo $accession does not start with GSE or equal E-MTAB-292...
    
  fi
  
done < $data/GSE124814_experiment_accessions.tsv

################################################################################
# GSE164677
################################################################################

if [[ -d $data/GSE164677 ]]; then

  echo Data for GSE164677 already exists and will not be downloaded.
  
else

  echo Downloading GSE164677
  mkdir $data/GSE164677
  GSE164677_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164677/suppl/GSE164677_Asian_MB_RNA-seq.txt.gz"
  curl -o $data/GSE164677/$(basename GSE164677_url) --silent $GSE164677_url
  sleep 5

fi

# check md5 sums of downloaded files
#echo Checking md5 sums of downloaded files ...
#md5sum --check --quiet data/md5_check_sums.tsv
