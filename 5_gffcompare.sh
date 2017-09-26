#!/usr/bin

## 2017-09-05
## RS Fraser, rfrase03@uoguelph.ca

## Compare the merged GTFs from our data to the reference dataset.

## gffcompare v0.10.1

## Variables
ref_gtf=/media/russ/data/cchf/genome/GRCh38.84.gtf

## Command
gffcompare \
	-r "$ref_gtf" \
    -G \
	-o compared \
    merged.gtf \
