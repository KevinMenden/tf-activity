#!/usr/bin/env python

## Script to intersect the region bed file of CAGE peaks
## with peak data from ChIP-seq experiments for confirmation
## of TF binding sites

import argparse
import pybedtools
import glob
import pandas as pd

## Some function
def get_tf_name(filename):
    sp = filename.split("/")
    tf_name = sp[len(sp)-1]
    tf_name = tf_name.split("_")
    return tf_name[0]

def bt_to_dataframe(bt):
    df = pd.read_table(bt.fn)
    return df

def get_bed_file_area(bf):
    """
    Calculate the total length of the region
    covered by a bed file
    :param bf: BedTool object
    :return: length of covered region
    """
    acc = 0
    df = bt_to_dataframe(bf)
    for i in range(len(df)):
        x = int(df.iloc[i, 2]) - int(df.iloc[i, 1])
        acc += x
    return acc


"""
Parser arguments
- region file containing the TSS peaks
- directory contaning peak files of TFs
"""
parser = argparse.ArgumentParser()
parser.add_argument("region_file", help="The region file of TSS peaks")
parser.add_argument("tf_dir", help="The directory containing the bed files of TF ChIP-seq peaks")
parser.add_argument("-out", help="The name of the output file")
args = parser.parse_args()
region_file = args.region_file
tf_dir = args.tf_dir
outfile = args.out

# Create BedTool object for region file
reg_bed = pybedtools.BedTool(region_file)
#reg_bed = reg_bed.sort()
#reg_bed = reg_bed.merge()

reg_bed_region = get_bed_file_area(reg_bed)
print("Region covered by TSS peaks: " + str(reg_bed_region))

# Get list of all TF bed files
tf_files = glob.glob(tf_dir + "/*.bed")

# Open output file
outfile = open(outfile, "w")
outfile.write("TF\tOverlaps\n")

tf_frames = []

# Loop over all TFs to find matches
for i,tff in enumerate(tf_files):
    tf_name = get_tf_name(tff)
    tf_bed = pybedtools.BedTool(tff)
    ovlp = reg_bed.intersect(tf_bed, u=True)
    ovlp_len = len(ovlp)
    df = pd.read_table(ovlp.fn, header=None, names=["chr", "start", "end", "strand"])
    df['TF'] = tf_name
    tf_frames.append(df)
    outfile.write(tf_name + "\t" + str(ovlp_len) + "\n")

outfile.close()

# Save TF-intersections
tfdf = pd.concat(tf_frames, ignore_index=True, axis=0)
csv_out_name = args.out
csv_out_name.replace(".txt", ".csv")
tfdf.to_csv(csv_out_name)