#!/usr/bin/sh

#RS Fraser 17-09-01

# Nanuq provides data in BAM format.
# HISAT2 requires FQ input; need to recreate the
# paired end files. GenomeQuebec recommends Picard
# tools for that task.

a=$1
base=${a%.bam}


java -Xmx20g -jar ~/java/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=$a FASTQ="$base".r1.fq SECOND_END_FASTQ="$base".r2.fq UNPAIRED_FASTQ="$base".unpaired.fq

## All the data would take up about 0.5Tb so gzip those suckahs

gzip "$base".r1.fq "$base".r2.fq "$base".unpaired.fq

