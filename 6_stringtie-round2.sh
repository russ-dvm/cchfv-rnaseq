#!/bin/bash

## RS Fraser, rfrase03@uoguelph.ca
## 2017-09-07

## Following the merge of GTF files generated by StringTie, re-run stringtie including the merged.gtf. This provides a uniform set of transcript assemblies for each lib.

## Stringtie v1.3.3b

## Variables

merged=/media/russ/data/cchf/5_merged/merged.gtf
output_dir=/media/russ/data/cchf/6_reestimate
output_deseq=/media/russ/data/cchf/6a_reestimate_for_deseq2
bam=$1
out="${bam%.bam}"
sample_dir="${out#HI*In*.}"

## Run on bam files.

echo "$bam for BALLGOWN"

stringtie \
    -e \
    -B \
    -p 8 \
    -G "$merged" \
    -o "$output_dir"/"$sample_dir"/"$out".gtf \
    $1

echo "$bam for DESEQ2"

stringtie \
    -e \
    -p 8 \
    -G "$merged" \
    -o "$output_deseq"/"$out".gtf \
    $1