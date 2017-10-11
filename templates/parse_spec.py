#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: parses dataframe to get productivity status                        #
#################################################################################

import math
import pandas as pd
import numpy as np
from collections import Counter

def _weighted_mu_std(values, weights):
    mu = np.average(values, weights=weights)
    sigma = np.average((values-mu)**2, weights=weights)
    return (mu, math.sqrt(sigma))

def _exclude_outliers(values, weights, m=4.0):
    values = np.array(values)
    weights = np.array(weights)
    avg, std = _weighted_mu_std(values, weights)
    selection = abs(values - avg) <= m * std
    return values[selection].tolist()


def _process_hdf(sample_name, region_file, refined_file, region, outdir=""):
    df = pd.read_hdf(region_file)
    check_prod = pd.read_hdf(refined_file)
    # 1-1 mapping
    assert (len(df) == len(check_prod))
    prods_only = df[(check_prod['v-jframe'] == 'In-frame') & (check_prod['stopcodon'] == 'No')]
    cdr3_lens = Counter(map(lambda x: len(x), prods_only[region].tolist()))

    sizes = cdr3_lens.keys()
    weights = map(lambda x : cdr3_lens[x], sizes)
    sizes = set(_exclude_outliers(sizes, weights))

    cdr3_lens_filtered = {k:v for k, v in cdr3_lens.iteritems() if k in sizes}

    with open(outdir+sample_name+"_{}_spec.csv".format(region), "w") as fp:
        header= ["length", "count"]
        fp.write(",".join(header))
        fp.write("\n")
        for k, v in cdr3_lens_filtered.items():
            fp.write("{},{}\n".format(k, v))



def process_samples(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone):
    fulldir1 = dir1+path_to_clone
    fulldir2 = dir2+path_to_clone
    fulldir3 = dir3+path_to_clone

    process_hdf(sample_name1, fulldir1+sample_name1+"_clones_seq.h5", fulldir1+sample_name1+"_refined_clones_annot.h5",
                'cdr3', dir1+"/diversity/")
    process_hdf(sample_name2, fulldir2+sample_name2+"_clones_seq.h5", fulldir2+sample_name2+"_refined_clones_annot.h5",
                'cdr3', dir2+"/diversity/")
    process_hdf(sample_name3, fulldir3+sample_name3+"_clones_seq.h5", fulldir3+sample_name3+"_refined_clones_annot.h5",
                'cdr3', dir3+"/diversity/")
