#!/usr/bin/env python

#================================================#
# Transform bed file into peak file with unique  #
# row ids                                        #
#================================================#

import sys

bed = sys.argv[1]
out = sys.argv[2]
bf = open(bed, "r")
of = open(out, "w")

bf_lines = bf.readlines()
bf.close()
for line in bf_lines:
    sp = line.split()
    row_id = "_".join(sp)
    of.write(row_id + "\t" + sp[0] + "\t" + sp[1] + "\t" + sp[2] + "\t" + sp[3] + "\n")
of.close()


