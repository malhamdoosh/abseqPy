#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: processes abundance file to prepare for a circos plot n R          #
#################################################################################

import pandas as pd
import re
from collections import defaultdict


def _get_canonical_family(name, family):
    # IGHV*/IGLV*
    if len(name) >= 5:
        # here, assume we can definitely find a V or J, otherwise .group() raises exception
        return family.search(name).group()
    return None


def _process_hdf(sample_name, fname, outdirname=""):
    # IGHVDJ or IGLVJ (Light chain has no D gene!)
    family = re.compile(r'IG(H[VDJ]|L[VJ])\d+')

    tally = defaultdict(lambda : defaultdict(int))

    df = pd.read_hdf(fname)

    for read in xrange(len(df)):
        v = _get_canonical_family(df.iloc[read]['vgene'], family)
        j = _get_canonical_family(df.iloc[read]['jgene'], family)
        if v is None or j is None:
            continue
        vgene = v.strip().upper()
        jgene = j.strip().upper()
        tally[vgene][jgene] += 1

    with open(outdirname + sample_name + "_vjassoc.csv", "w") as fp:
        header = ["from", "to", "value"]
        fp.write(",".join(header))
        fp.write("\n")

        for vgene, dic in tally.items():
            for jgene, value in dic.items():
                fp.write("{},{},{}\n".format(vgene, jgene, value))



def process_samples(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone):
    fulldir1 = dir1+path_to_clone
    fulldir2 = dir2+path_to_clone
    fulldir3 = dir3+path_to_clone

    process_hdf(sample_name1, fulldir1+sample_name1+"_refined_clones_annot.h5", dir1+"/abundance/")
    process_hdf(sample_name2, fulldir2+sample_name2+"_refined_clones_annot.h5", dir2+"/abundance/")
    process_hdf(sample_name3, fulldir3+sample_name3+"_refined_clones_annot.h5", dir3+"/abundance/")
