#!/bin/bash

base_dir="/root/microbiome_meta/project5"
gsa_base="https://download.cncb.ac.cn/gsa/CRA003252"

for type in 16S ITS
do
  sample_sheet="${base_dir}/${type}_samples.tsv"
  [ -f "$sample_sheet" ] || continue

  tail -n +2 "$sample_sheet" | while IFS=$'\t' read accession group
  do
    # Make sure group name is lowercase (control/treated)
    group=$(echo "$group" | tr '[:upper:]' '[:lower:]')
    out_dir="${base_dir}/${type}/${group}"
    mkdir -p "$out_dir"

    echo "🔽 Downloading ${accession}_f1.fastq.gz to $out_dir/${accession}_1.fastq.gz"
    wget -c "$gsa_base/$accession/${accession}_f1.fastq.gz" -O "$out_dir/${accession}_1.fastq.gz"

    echo "🔽 Downloading ${accession}_r2.fastq.gz to $out_dir/${accession}_2.fastq.gz"
    wget -c "$gsa_base/$accession/${accession}_r2.fastq.gz" -O "$out_dir/${accession}_2.fastq.gz"
  done
done
