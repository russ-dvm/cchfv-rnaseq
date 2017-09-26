#!/usr/bin/python

from __future__ import print_function
import argparse
import os


parser = argparse.ArgumentParser(
    description="Create a summary from the summary files of Hisat2 (not the new-summary) or STAR.",
    epilog="RS Fraser, 17-09-07"
)
parser.add_argument(
    "summary",
    help="Summary file from HISAT2 or STAR"
)
parser.add_argument(
    "--hisat",
    help="File comes from HISAT2",
    action="store_true"
)
parser.add_argument(
    "--star",
    help="File comes from STAR",
    action="store_true"
)

args = parser.parse_args()
in_file = open(args.summary, 'r')


def hisat(file):
    print("lib", "aligner", "category", "value", sep="\t")
    a, base, extension = file.name.split(".")
    lines = file.readlines()
    total = lines[0].split()[0]
    print("hisat", base, "total", total, sep="\t")
    paired = lines[1].split()[0]
    print("hisat", base, "paired", paired, sep="\t")
    paired_concord_0 = lines[2].split()[0]
    print("hisat", base, "paired_concord_0", paired_concord_0, sep="\t")
    paired_concord_1 = lines[3].split()[0]
    print("hisat", base, "paired_concord_1", paired_concord_1, sep="\t")
    paired_concord_2 = lines[4].split()[0]
    print("hisat", base, "paired_concord_2", paired_concord_2, sep="\t")
    discord_1 = lines[7].split()[0]
    print("hisat", base, "discord_1", discord_1, sep="\t")
    rate = lines[14].split()[0]
    rate = rate[:-1]
    print("hisat", base, "rate", rate, sep="\t")


def star(file):
    pass


if args.hisat:
    hisat(in_file)
elif args.star:
    star(in_file)
else:
    print("Please specify the aligner used (--hisat or --star)")
