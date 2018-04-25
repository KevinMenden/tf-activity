#!/usr/bin/env python

"""
Filter motif instance results to create list of TF targets
"""

# Imports

import argparse
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests

def append_annotations(df, anno):
    """
    Append the annotation column from the anno df to the target df
    :param df: the target-gene mapping df
    :param anno: df containing the annotation
    :return: target-gene mapping df with annotation
    """
    targets = df['Target']
    anno_list = []
    for tar in targets:
        tmp = anno.loc[anno['PeakID'] == tar]
        anno_list.append(list(tmp['Annotation']))
    anno_list_flat = [y for x in anno_list for y in x]
    df['Annotation'] = anno_list_flat
    return df


def calculuate_hits(enriched_tfs, inst_df):
    """
    Calculate number of gene hits for each TF and gene pair
    :param enriched_tfs: DF of enrichment results
    :param inst_df: DF of motif instance finding results
    :return: DataFrame containing Target, TF and Number of hits
    """
    df_dict = {}
    for tf in enriched_tfs:
        hits = inst_df.loc[inst_df['motif'] == tf]
        targets = pd.Series(hits['peak'])
        # Count number of hits for each target
        utargets = set(targets)
        tlist = list(targets)
        ut_dict = {}
        for ut in utargets:
            no_hits = tlist.count(ut)
            ut_dict[ut] = no_hits

        # Create data frame containing target, hits and TF
        ut_df = pd.DataFrame.from_dict(ut_dict, orient='index')
        ut_df['TF'] = pd.Series([tf] * len(ut_df), index=ut_df.index)
        ut_df.index.name = "Target"
        ut_df.reset_index(inplace=True)

        # Add to dict
        df_dict[tf] = ut_df
    mdata = pd.DataFrame(pd.concat(df_dict))
    return mdata


# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("enrichment", help="Result from homer2 enrichment analysis")
parser.add_argument("inst", help="Result from home2 find motifs analysis")
parser.add_argument("anno_peaks", help="The annotate CAGE peaks")
parser.add_argument("--pval_cutoff", help="The adj. P-value cutoff. Default 0.05", default=0.05)
parser.add_argument("--log_odds", help="The log-odds cutoff for hit decision. Defaul 9.0", default=9.0)
args = parser.parse_args()
enr = args.enrichment
inst = args.inst
anno = args.anno_peaks

# Load instance result dataframe
inst_df = pd.read_table(inst, header=None)
inst_df.columns = ['peak', 'location', 'instance', 'motif', 'strand', 'log_odds']
# Load enrichment dataframe
enr_df = pd.read_table(enr)
# Load annotation dataframe
anno_df = pd.read_table(anno)
anno_df = anno_df.rename(columns={anno_df.columns.tolist()[0]: 'PeakID'})
anno_df = anno_df[['PeakID', 'Annotation']]

# Set Cutoffs
pval_cutoff = args.pval_cutoff
log_odds_cutoff = args.log_odds

# Adjust p-values and get significant enrichment results
pvals = enr_df['p-value']
padj = multipletests(pvals, method="fdr_bh")[1]
enr_df['padj'] = pd.Series(padj)
enr_df_sig = enr_df.loc[enr_df['padj'] <= pval_cutoff]
enriched_tfs = enr_df_sig['Motif Name']
# Get above-threshod motif instances
inst_df_sig = inst_df.loc[inst_df['log_odds'] >= log_odds_cutoff]

## Calculate number of hits per TF and target and save
mdata = calculuate_hits(enriched_tfs, inst_df_sig)
mdata.to_csv("tf_target_mapping.txt", sep="\t")

# Merge with annotation and save
target_anno_df = append_annotations(mdata, anno_df)
target_anno_df.to_csv("tf_target_mapping_annot.txt", sep="\t")


