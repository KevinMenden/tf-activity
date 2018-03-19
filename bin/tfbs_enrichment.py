# Transcription Factor Binding Site enrichment analysis on CAGE peaks

# Imports
import os
import argparse
import datetime
import pandas as pd
from Bio import SeqIO
from Bio import motifs


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


def extract_extended_regions(tss_coords, extension_range, genome_dict):
    """
    Extract sequence of extended regions around TSS peaks from CAGE
    :param tss_coords: the coordinates in format chr_start_end_strand
    :param extension_range: the range of extension in bp
    :param genome_dict: the dictionary containing the genome sequences
    :return: Series of extended regions as Seq objects
    """
    tss_regions = []
    for tssc in tss_coords:
        tss = tssc.split("_")
        chr = tss[0]
        start = int(tss[1]) - (2*extension_range) - 1
        end = int(tss[2]) + extension_range - 1
        strand = tss[3]
        rec = genome_dict[chr]
        tss_reg = rec.seq[start:end]
        tss_regions.append(tss_reg)
    return pd.Series(tss_regions)

def count_gen(iter):
    """
    Count the items in the iterator object with the second tuple element
    greate than 0
    :param iter:
    :return:
    """
    count = 0
    for elem in iter:
        if elem[1] > 0:
            count +=1
    return count


def set_alphabet_to_motif_alphabet(seqs, mot):
    """
    Change alphabet of Seq objects to alphabet of motif object
    :param seqs: Series of Seq objects
    :param mot: a motif object
    :return: updated Series of Seq objects
    """
    for seq in seqs:
        seq.alphabet = mot.alphabet
    return seqs


def write_config_file(args, analysis_dir, ts):
    """
    Write config file
    :param args:
    :param analysis_dir:
    :return:
    """
    nf = open(analysis_dir + "/config.cfg", 'w')
    nf.write("Time:\t\t" + ts + "\n")
    nf.write("Peak file\t\t" + args.tss_peaks + "\n")
    nf.write("Genome file\t\t" + args.genome + "\n")
    nf.write("Range\t\t" + str(args.range) + "\n")
    nf.write("Motif file\t\t" + args.motifs + "\n")
    nf.write("Output directory\t\t" + args.output_dir + "\n")
    nf.write("Log Odds Score\t\t" + str(args.log_odds_score) + "\n")
    nf.write("P-value\t\t" + str(args.pval) + "\n")
    nf.write("LFC cutoff\t\t" + str(args.lfc_cutoff) + "\n")
    nf.write("TFBS searched\t\t" + str(args.find_instances) + "\n")
    nf.close()

def find_binding_sites(mots, seqs, tss_coords ,log_odds_score, analysis_dir):
    print("Calculating odds scores ...")
    tf_dict = {}
    background = {'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}
    pseudocounts = {'A': 0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
    no_motifs = len(mots)
    for i,m in enumerate(mots):
        score_acc = 0
        pwm = m.counts.normalize(pseudocounts=pseudocounts)
        pssm = pwm.log_odds(background)
        print("Processing " + m.name + "\t" + str(i+1) + " out of " + str(no_motifs))
        nf = open(analysis_dir + "/" + m.name + "_hits.txt", 'w')
        nf.write("TF\tTSS_Peak\tHits\n")
        for idx, seq in enumerate(seqs):
            hits = count_gen(pssm.search(seq, log_odds_score))
            score_acc += hits
            nf.write(m.name + "\t" + tss_coords[idx] + "\t" + str(hits) + "\n")
        tf_dict[m.name] = score_acc
        nf.close()

    return tf_dict

#####################
### MAIN FUNCTION ###
#####################
if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("tss_peaks", help="The file containing the DE table for the TSS peaks")
    parser.add_argument("-g", "--genome", help="The genome in fasta format",
                        default="/home/kevin/resources/genomes/GRCh38_v27_gencode/GRCh38.primary_assembly.genome.fa")
    parser.add_argument("-r", "--range", help="The range used to extend the peaks up- and downstream (twice the length). Default = 300",
                        default=300)
    parser.add_argument("-m", "--motifs", help="The file containing the motifs of interest in jaspar format",
                        default="/home/kevin/resources/TF_interactions/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
    parser.add_argument("-o", "--output_dir", help="The directory to store the results in. Default: working directory.",
                        default=".")
    parser.add_argument("-lod", "--log_odds_score", help="The log odd score cutoff. Default = 10.0",
                        default=10.0)
    parser.add_argument("-pval", help="P-value cutoff for TSS peaks. Default = 0.05",
                        default=0.05)
    parser.add_argument("-lfc", "--lfc_cutoff", help="Log-fold change cutoff for TSS peaks. Default = 1",
                        default=1)
    parser.add_argument("--find_instances", help="Look for instances of the motif in the extracted sequence. Default is false. Can be slow.",
                        action="store_true", default=False)
    parser.add_argument("-p", "--project", help="Optional name to add to the directory and file names", default="")

    args = parser.parse_args()
    peak_file = args.tss_peaks
    genome_file = args.genome
    extension_range = int(args.range)
    jasp_motifs = args.motifs
    output_dir = args.output_dir
    log_odds_score = float(args.log_odds_score)
    pval = float(args.pval)
    lfc_cutoff = float(args.lfc_cutoff)
    project = str(args.project)
    if not project == "":
        project = project + "_"  # add this for better readability of project name given

    # Create new directory
    ts = datetime.datetime.now().isoformat().split(".")[0].replace(":", "")
    analysis_dir = output_dir + "/tfbs_enrichment_analysis_" + project + ts
    os.mkdir(analysis_dir)

    # Create config file
    write_config_file(args, analysis_dir, ts)


    fasta_filename = analysis_dir + "/extended_peaks_seqs_" + project + ts + ".fasta"

    # Load the genome
    print("\n#==== Loading genome ====# \n" + genome_file)
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    print("Genome loaded\n")

    # Read peak file peak file
    peak_df = pd.read_table(peak_file)
    peak_df = peak_df.rename(columns={peak_df.columns.tolist()[0]: 'tss_coordinates'})

    # Only consider significant, up-regulated peaks
    peak_df = peak_df.loc[peak_df['padj'] <= pval]
    peak_df = peak_df.loc[peak_df['log2FoldChange'] >= lfc_cutoff]
    peak_df.index = pd.RangeIndex(peak_df.shape[0])

    tss_coords = peak_df['tss_coordinates']

    # Extract sequence of extended regions
    print("Extracting regions ...")
    tss_regions = extract_extended_regions(tss_coords, extension_range, genome_dict)
    tss_coords = pd.Series(tss_coords)

    # Save region sequences as fasta file
    print("Writing extended sequences to fasta ...")
    write_sequence_fasta(fasta_filename, tss_regions, tss_coords)

    if args.find_instances:
        # Load motifs
        print("Loading motifs ...")
        fh = open(jasp_motifs)
        tf_motifs = motifs.parse(fh, 'jaspar')

        # Update alphabet of sequences
        tss_regions = set_alphabet_to_motif_alphabet(tss_regions, tf_motifs[0])

        print("Searching for binding sites")
        tf_dict = find_binding_sites(tf_motifs, tss_regions, tss_coords, log_odds_score, analysis_dir)
        tf_dict_df = pd.Series(tf_dict).to_frame()
        tf_dict_df.columns = ['Score']
        tf_dict_df.to_csv(analysis_dir+"/tf_scores_"+ts+".csv", sep="\t")
