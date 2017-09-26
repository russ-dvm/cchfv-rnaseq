#!/usr/bin

## 2017-09-05
## RS Fraser, rfrase03@uoguelph.ca

## Assembles and quantifies the expressed genes and transcripts based on the output of HISAT2.

## Stringtie v1.3.3b

## Template file name: HI.4077.002.Index_12.Lib_18.bam

## Variable definitions
bam=$1
gtf=/media/russ/data/cchf/genome/GRCh38.84.gtf
output_dir=/media/russ/data/cchf/4_assembled
out="${bam%.bam}"
lib="${out#HI*In*.}"

echo $out, $lib

## Command. Remaining options are left at default values. 
## -c: Sets the minimum read coverage allowed for the predicted transcripts. A transcript with a lower coverage than this value is not shown in the output. Default: 2.5
## -l: Sets <label> as the prefix for the name of the output transcripts. Default: STRG

stringtie \
    -p 8 \
	-G "$gtf" \
	-o "$output_dir"/"$out".gtf \
	-l "$lib" \
    --rf \
    -c 2.5 \
	-A "$output_dir"/"$out".gene_abund.tab \
	-C "$output_dir"/"$out".cov_refs.tab \
    $bam

