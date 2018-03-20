#!/usr/bin/env python


## Use an "extended_peaks_seq" file to create a bed file defining the regions
##

import argparse
import os

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("peak_fasta", help="The extended peaks fasta file.")
parser.add_argument("out", help="Name of the output bedfile")
parser.add_argument("-range", help="The range for extension. Default = 300", default=300)
args = parser.parse_args()

peak_fasta = args.peak_fasta
outfile = args.out
rng = int(args.range)

outfile = open(outfile, "w")
fasta = open(peak_fasta, "r")

fasta_lines = fasta.readlines()
fasta.close()
peaks = [x for x in fasta_lines if x.startswith(">")]

for peak in peaks:
    peak_range = peak.split("_")
    chr = peak_range[0].replace(">", "")
    start = str(int(peak_range[1]) - (2*rng))
    end = str(int(peak_range[2]) + rng)
    strand = peak_range[3]
    outfile.write(chr + "\t" + start + "\t" + end + "\t" + strand)

outfile.close()