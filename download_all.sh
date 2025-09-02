#!/bin/bash

base_dir="/root/microbiome_meta"

for project in project1 project2 project3 project4
do
  for type in 16S ITS
  do
    sample_sheet="${base_dir}/${project}/${type}_samples.tsv"
    # Set output directory for this project/type
    out_dir="${base_dir}/${project}/${type}"
    mkdir -p "$out_dir"

    # Skip if sample sheet doesn't exist
    [ -f "$sample_sheet" ] || continue

    # Loop through Run_accession (skip header)
    tail -n +2 "$sample_sheet" | cut -f1 | while read accession
    do
      echo "Downloading $accession for $project $type"
      fasterq-dump "$accession" -O "$out_dir" --split-files
    done
  done
done
