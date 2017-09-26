#!/bin/bash
#RSF, rfrase03@uoguelph.ca

#run fastqc on all the fastq files.

raw=/home/russ/cchf/1_raw_data

for file in "$raw"/*r1.fq.gz
do
	read1="${file#/home/russ/cchf/1_raw_data/}"
	group_num="${read1%1.fq.gz}"
	read2="$group_num"2.fq.gz
	fastqc --nogroup --noextract "$raw"/"$read1" "$raw"/"$read2" -o /home/russ/cchf/2_fastqc/
done
