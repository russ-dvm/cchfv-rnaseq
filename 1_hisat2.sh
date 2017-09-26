#!/usr/bash

## 2017-09-04
## RS Fraser, rfrase03@uoguelph.ca

## Command used to align f/r reads to human genome.

## HISAT2 version 2.1.0


## sample filename: HI.4077.001.Index_2.Lib_1.r1.fq.gz


f=$1
base="${f%1.fq.gz}"
r="$base"2.fq.gz
lib="${base%.r}"
lib="${lib#HI*x_*.}"
output="${base%.r}"

index_dir=/home/russ/cchf/genome/grch38_snp_tran
index=genome_snp_tran


output_dir_aligned=/home/russ/cchf/3_aligned
output_dir_unaligned=/home/russ/cchf/3b_unaligned


## --known-splicesite-infile: not needed; included in the hisat2 index
## --transcriptome-mapping-only - worth considering? Do we want novel transcripts?

hisat2 \
	-p 8 \
	--dta \
	-t \
	--rg-id "$lib" \
	--rg PL:ILLUMINA \
	--un-conc-gz "$output_dir_unaligned"/"${base%.r}".unaligned.fq.gz \
	--summary-file "$output_dir_aligned"/summary."$lib".txt \
	--rna-strandness RF \
	-x "$index_dir"/"$index" \
	-1 $f \
	-2 $r \
	-S "$output_dir_aligned"/"$output".sam \


## Clean up the mega huge sam files. Could probably pipe the output of hisat into samtools but for some reason doesn't seem to be a thing in the documentation.
samtools sort -@ 8 -o "$output_dir_aligned"/"$output".bam "$output_dir_aligned"/"$output".sam
rm "$output_dir_aligned"/"$output".sam