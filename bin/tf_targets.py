#!/usr/bin/env python

"""
Filter motif instance results to create list of TF targets
"""

# Imports

import argparse
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests

############ END FUNCTIONS #################


# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("enrichment", help="Result from homer2 enrichment analysis")
parser.add_argument("inst", help="Result from home2 find motifs analysis")
args = parser.parse_args()
enr = args.enrichment
inst = args.inst

inst_df = pd.read_table(inst, header=None)
inst_df.columns = ['peak', 'location', 'instance', 'motif', 'strand', 'log_odds']

enr_df = pd.read_table(enr)

# Cutoffs
pval_cutoff = 0.001
log_odds_cutoff = 9.0

# Adjust p-values and get significant enrichment results
pvals = enr_df['p-value']
padj = multipletests(pvals)[1]
enr_df['padj'] = pd.Series(padj)
enr_df_sig = enr_df.loc[enr_df['padj'] <= pval_cutoff]
enriched_tfs = enr_df_sig['Motif Name']

# Get above-threshod motif instances
inst_df_sig = inst_df.loc[inst_df['log_odds'] >= log_odds_cutoff]

######
## Calculate number of hits per TF and target
######
df_dict = {}
for tf in enriched_tfs:
    hits = inst_df_sig.loc[inst_df_sig['motif'] == tf]
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

# Concatenate everything
mdata = pd.DataFrame(pd.concat(df_dict))
mdata.to_csv("tf_target_mapping.txt")

######
## Convert DF to BED format for annotation
######


