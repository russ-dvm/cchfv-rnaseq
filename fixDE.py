#!/usr/bin/python


from __future__ import print_function
import argparse
from pyensembl import EnsemblRelease

parser = argparse.ArgumentParser(description = "This script is designed to ensure that each MSTRG has an ENS gene and transcript id, and gene name. Meant to be run on the output of Stringtie (ie a gtf file). NOTE: The file must be SORTED before running (sort -k 8)", epilog = "RS Fraser, 2017-10-9")

parser.add_argument("gtf", help = "A GTF file created following the Stringtie 'merge' step.")


args = parser.parse_args()

gtf = open(args.gtf, 'r')

data = EnsemblRelease(90)

## Gonna make a big ol' list of gene ids so that we don't have
## multiple entries for the same MSTRG id.
g_id = []
mstrg_list = []

for line in gtf:
    if "#" in line:
        continue
    else:
        line = line.rstrip()
        line_fields = line.split("\t")
        tran = line_fields[2]

        ## only consider transcript entries (avoid exons)
        if "transcript" in tran:
            t_id = line_fields[8].split(";")[1]
            g_name = line_fields[8].split(";")[2]
            g_id_fresh = line_fields[8].split(";")[0]

            ## Only print one gene for multiple transcripts - once we've printed the info for a gene, add that gene id to the gene_id list.
            if g_id_fresh in g_id:
                continue 

            ## But, we need to check for an ENS id in all transcripts for the same gene - occasionally the first hit will be a novel transcript
            ## Check if the transcript is known, if not, add it to a tracking list
            if "MSTRG" in t_id:              
                if g_id_fresh not in mstrg_list and len(mstrg_list) != 0:
                    mstrg = line2.split("\t")[8].split(";")[0].split("\"")[1]
                    mstrg_t_id = line2.split("\t")[8].split(";")[1].split("\"")[1]
                    print(mstrg, mstrg, mstrg_t_id, mstrg, sep = "\t")
                    mstrg_list = []
                    mstrg_list.append(g_id_fresh)
                
                else:
                    mstrg_list.append(g_id_fresh)


            ## Check if our gene id is already in the MSTRG tracker
            if g_id_fresh in mstrg_list:

                ## If subsequent entries for the MSTRG gene have an ensembl transcript, print out the correct info
                if "ENS" in t_id:
                    mstrg = g_id_fresh.split("\"")[1] 
                    g_name = g_name.split("\"")[1]
                    t_id = t_id.split("\"")[1]

                    ensembl = data.transcript_by_id(t_id)
                    ensembl_id = ensembl.gene_id
                    ensembl_name = ensembl.transcript_name                
                    print(mstrg, ensembl_id, t_id, g_name, sep = "\t")

                    ##  kick future genes with this id out of the loop
                    g_id.append(g_id_fresh)
                    ## and squash the mstrg list
                    mstrg_list = []

                ## if no ENS id is found, save the info for the next round
                elif "MSTRG" in t_id:
                    line2 = line

            # if the id hasn't been found yet, and is NOT in the mstrg_list, then:
            else:
                ## print off whatever was left over in the mstrg list
                if len(mstrg_list) != 0:
                    mstrg = line2.split("\t")[8].split(";")[0].split("\"")[1]
                    mstrg_t_id = line2.split("\t")[8].split(";")[1].split("\"")[1]
                    print(mstrg, mstrg, mstrg_t_id, mstrg, sep = "\t")
                    mstrg_list = []

                ## then, print off the info for this gene
                mstrg = g_id_fresh.split("\"")[1] 
                g_name = g_name.split("\"")[1]
                t_id = t_id.split("\"")[1]

                ensembl = data.transcript_by_id(t_id)
                ensembl_id = ensembl.gene_id
                ensembl_name = ensembl.transcript_name                    
                print(mstrg, ensembl_id, t_id, g_name, sep = "\t")
                g_id.append(g_id_fresh)
