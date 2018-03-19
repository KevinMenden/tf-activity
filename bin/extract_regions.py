#!/usr/bin/env python

# Transcription Factor Binding Site enrichment analysis on CAGE peaks

# Imports
import argparse
import datetime
import pandas as pd
from Bio import SeqIO


###### FUNCTION SECTION ######

# Create DataFrame
def write_sequence_fasta(fasta_file, tss_regions_sequences, tss_names):
    """
    Write fasta format file with the extendede sequences from the TSS peaks
    :param fasta_file: name of the output file
    :param tss_regions_sequences: the extended TSS region sequences
    :param tss_names: the names of the TSS peaks (or coordinates)
    :return: null
    """
    nf = open(fasta_file, 'w')
    for i in range(len(tss_regions_sequences)):
        tmp_name = tss_names[i]
        tmp_reg = str(tss_regions_sequences[i])
        nf.write(">" + tmp_name + "\n")
        nf.write(tmp_reg + "\n")
    nf.close()


def extract_extended_regions(peak_file, extension_range, genome_dict):
    """
    Extract sequence of extended regions around TSS peaks from CAGE
    :param peak_file: the bed file containing the peaks
    :param extension_range: the range of extension in bp
    :param genome_dict: the dictionary containing the genome sequences
    :return: Series of extended regions as Seq objects
    """
    peaks = open(peak_file, "r").readlines()
    tss_regions = []
    tss_coords = []
    for peak in peaks:
        tss_list = peak.split("\t")
        chr = tss_list[0]
        start = int(tss_list[1]) - (2*extension_range) - 1
        end = int(tss_list[2]) + extension_range - 1
        strand = tss_list[3].rstrip()
        tss_coords.append("_".join(tss_list).rstrip())
        rec = genome_dict[chr]
        if strand == '+':
            tss_reg = rec.seq[start:end]
        else:
            tss_reg = rec.seq[start:end].reverse_complement()
        tss_regions.append(tss_reg)
    return [pd.Series(tss_regions), pd.Series(tss_coords)]


#####################
### MAIN FUNCTION ###
#####################
if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("tss_peaks", help="The bed file containing the DE table for the TSS peaks")
    parser.add_argument("-g", "--genome", help="The genome in fasta format",
                        default="/home/kevin/resources/genomes/GRCh38_v27_gencode/GRCh38.primary_assembly.genome.fa")
    parser.add_argument("-r", "--range", help="The range used to extend the peaks up- and downstream (twice the length). Default = 300",
                        default=300)
    parser.add_argument("-p", "--project", help="Optional name to add to the directory and file names", default="")

    args = parser.parse_args()
    peak_file = args.tss_peaks
    genome_file = args.genome
    extension_range = int(args.range)
    project = str(args.project)
    if not project == "":
        project = project + "_"  # add this for better readability of project name given

    fasta_filename = "extended_peaks_seqs_" + project + ".fasta"

    # Load the genome
    print("\n#==== Loading genome ====# \n" + genome_file)
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    print("#==== Genome loaded ====# \n")


    # Extract sequence of extended regions
    print("Extracting regions ...")
    [tss_regions, tss_coords] = extract_extended_regions(peak_file, extension_range, genome_dict)


    # Save region sequences as fasta file
    print("Writing extended sequences to fasta ...")
    write_sequence_fasta(fasta_filename, tss_regions, tss_coords)

