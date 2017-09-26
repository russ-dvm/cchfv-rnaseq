#!/usr/bin

## 2017-09-05
## RS Fraser, rfrase03@uoguelph.ca

## Merges the output of all stringtie assemblies into a single GTF.

## Stringtie v1.3.3b

## Variables

gtf=/media/russ/data/cchf/genome/GRCh38.84.gtf
output_dir=/media/russ/data/cchf/5_merged


## Create file list
echo "Creating list of gtfs..."
for x in *gtf
do
    echo $x
	readlink -f $x >> gtf.list
done

## Command
echo "Stringtie-ing"
stringtie \
	--merge \
	-p 8 \
	-G "$gtf" \
	-o "$output_dir"/merged.gtf \
	gtf.list