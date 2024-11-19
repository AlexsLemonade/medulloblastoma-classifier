#!/usr/bin/bash

set -euo pipefail

file=$1
output_dir=$2

mkdir -p ${output_dir}
gunzip -c ${file} > ${output_dir}/temp.csv

previous_end=1

# count each column prefix to define number of columns for each sample ID
head -1 ${output_dir}/temp.csv | tr ',' '\n' | sed '1d' | cut -f1 -d'_' | uniq -c | sed 's/^\s\+//' | while read n_cols sample_id; do

  # update the start and end column based on number of columns for this sample
  this_start=$((${previous_end} + 1))
  this_end=$((${this_start} + ${n_cols} - 1))

  # cut the gene column (1) and this sample's columns and write to output
  gunzip -c ${file} | cut -f1,${this_start}-${this_end} -d',' | tr ',' '\t' > ${output_dir}/${sample_id}.tsv

  # update previous_end for next sample
  previous_end=${this_end}

done

rm -f ${output_dir}/temp.csv
