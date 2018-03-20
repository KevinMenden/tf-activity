#!/usr/bin/env python

"""
Find instances of TF motifs in the provided sequences
"""

# Imports
import argparse
from Bio import SeqIO
from Bio import motifs


########## FUNCTION SECTION ###############
def set_alphabet_to_motif_alphabet(seq_dict, mot):
    """
    Change alphabet of Seq objects to alphabet of motif object
    :param seqs: Series of Seq objects
    :param mot: a motif object
    :return: updated Series of Seq objects
    """
    for elem in seq_dict:
        seq_rec = seq_dict[elem]
        seq = seq_rec.seq
        seq.alphabet = mot.alphabet
        seq_dict[elem] = seq
    return seq_dict


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
            count += 1
    return count


def find_binding_sites(mts, seqs, log_odds_score=8.0):
    """
    Find instances of motifs in sequences
    :param mts:
    :param seqs:
    :param log_odds_score:
    :return:
    """
    print("Calculating odds scores ...")
    background = {'A': 0.3, 'C': 0.2, 'G': 0.2, 'T': 0.3}
    pseudocounts = {'A': 0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
    no_motifs = len(mts)
    for i, m in enumerate(mts):
        pwm = m.counts.normalize(pseudocounts=pseudocounts)
        pssm = pwm.log_odds(background)
        print("Processing " + m.name + "\t" + str(i + 1) + " out of " + str(no_motifs))
        nf = open(m.name + "_hits.txt", 'w')
        nf.write("TF\tSequence\tHits\n")
        for elem in seqs:
            seq = seqs[elem]
            hits = count_gen(pssm.search(seq, log_odds_score))
            nf.write(m.name + "\t" + elem + "\t" + str(hits) + "\n")
        nf.close()


############ END FUNCTIONS #################


# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("fasta", help="The fasta sequence file in which to search for motif instances")
parser.add_argument("motifs", help="The file containing the TF motifs in .jaspar motifs")
args = parser.parse_args()
fasta = args.fasta
mots = args.motifs

############# Find instances ###############
# Parse motifs
fh = open(mots)
tf_motifs = motifs.parse(fh, 'jaspar')
# TODO delete short-listing
tf_motifs = tf_motifs[0:5]
# Parse target sequences
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
seq_dict = set_alphabet_to_motif_alphabet(seq_dict, tf_motifs[0])

# Find motif instances
find_binding_sites(tf_motifs, seq_dict)


############## EOF ##########################