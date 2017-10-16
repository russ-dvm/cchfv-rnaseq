#!/usr/bin/python

from __future__ import print_function
import sys

fpkm = open(sys.argv[1], 'r')

print("library", "gene_id", "transcript", "FPKM", "TPM", sep = "\t")

lib = sys.argv[1]
start_in = lib.index("Lib_")
end_in = lib.index(".trimm")
lib = lib[start_in:end_in]


for line in fpkm:
    if "#" in line:
        continue
    line = line.rstrip()
    lineFields = line.split("\t")
    if "transcript" in lineFields[2]:
        info = lineFields[8].split(";")
        for field in info:
            if field.startswith("gene_id"):
                gene = field.split("\"")[1]
            elif field.startswith(" transcript_id"):
                trans = field.split("\"")[1]
            elif field.startswith(" FPKM"):
                fpkm_count = field.split("\"")[1]
            elif field.startswith(" TPM"):
                tpm = field.split("\"")[1]

        print(lib, gene, trans, fpkm_count, tpm, sep = "\t")
