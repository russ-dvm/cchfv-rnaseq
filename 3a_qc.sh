#!/usr/bash

## RS Fraser, rfrase03@uoguelph.ca
## 2017-09-05

## Set of commands to run basic QC on aligned BAMs prior to proceeding.
## Need to be in the directory with .bam files.

## Shared output directory
output_dir=/media/russ/data/cchf/3c_qc

## Check 5' - 3' bias - takes a million years to run through all bam files.
geneBody_coverage.py \
	-r /media/russ/data/cchf/ref_files/hg38.HouseKeepingGenes.bed \
	-i /media/russ/data/cchf/3_aligned/ \
	-o "$output_dir"/genebody/genebody


## Check splice junction detection (a way of determining whether enough sequencing depth was used to detect a) known splice junctions and/or b) novel splice junctions). Look to see whether the lines are approaching a plateau to determine whether enough depth was achieved.


for x in *bam
do

	out="${x%.bam}"
	out="${out#H*x_*.}"

	junction_saturation.py \
		-r /media/russ/data/cchf/ref_files/hg38.GENCODE26.bed \
		-i $x \
		-o "$output_dir"/jnc_sat/"$out"


	## Alignment size - does the average alignment size match the size selection from library prep?
	inner_distance.py \
		-i $x \
		-r /media/russ/data/cchf/ref_files/hg38.GENCODE26.bed \
		-o "$output_dir"/align_size/"$out"

done