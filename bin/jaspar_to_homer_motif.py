#!/usr/bin/env python

# Transform motifs in JASPAR format to motifs in HOMER format

import argparse
from Bio import motifs

# Function to normalize nucleotide counts
def norm_motif(line):
    fl = [float(i) for i in line]
    ls = sum(fl)
    nl = [i/ls for i in fl]
    nl = [str(i) for i in nl]
    return nl


# Function to create Homer motif from jaspar motif
def make_homer_motif(mot, cutoff):
    header = ">" + str(mot.consensus) + "\t" + mot.name + "\t" + cutoff + "\n"
    # Transfer to transfac format
    mot = mot.format("transfac")
    new_lines = ""
    mot_lines = mot.splitlines()
    for i in range(1,len(mot_lines)-2):
        line = mot_lines[i]
        sp = line.split()
        nl = norm_motif(sp[1:5])
        nl = "\t".join(nl)
        new_lines = new_lines + nl + "\n"
    new_mot = header + new_lines
    return new_mot

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("jaspar_motifs", help="The file containing the DE table for the TSS peaks")
parser.add_argument("-out", help="Optional name of output file", default="homer_motifs.homer")
parser.add_argument("-c", help="Log odds cutoff. Defaults to 9.0", default="9.0")
args = parser.parse_args()
jasp_motifs = args.jaspar_motifs
out_name = args.out
lod = str(args.c)

# Read jaspar motifs
fh = open(jasp_motifs)
jaspar = motifs.parse(fh, 'jaspar')

# Transform them and write them to new file
nf = open(out_name, "w")
for m in jaspar:
    homer_mot = make_homer_motif(m,lod)
    nf.write(homer_mot)

nf.close()
