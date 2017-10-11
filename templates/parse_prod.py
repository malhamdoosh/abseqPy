#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: parses dataframe to get productivity status                        #
#################################################################################

import pandas as pd
from collections import Counter

def _process_hdf(sample_name, fname, outdir=""):
    df = pd.read_hdf(fname)
    vjframe = Counter(df['v-jframe'].tolist())
    stopcod = Counter(df['stopcodon'].tolist())
    stopcod_inframe = len(df[(df['v-jframe'] == 'In-frame') & (df['stopcodon'] == 'Yes')])
    outframe_nostop = len(df[(df['v-jframe'] == 'Out-of-frame') & (df['stopcodon'] == 'No')])
    both = len(df[(df['v-jframe'] == 'Out-of-frame') & (df['stopcodon'] == 'Yes')])
    prod_reads = len(df[(df['v-jframe'] == 'In-frame') & (df['stopcodon'] == 'No')])

    # percentage calculated from total of Out-Of-Frame + In-frame
    frame_sample_total = sum(vjframe.values())
    # percentage calculated from total of stopcodon + no stop codon
    cod_sample_total = sum(stopcod.values())
    # percentage calculated from total sample size
    total_size = len(df)

    res = {
        'Productivity': ["Productive", "Unproductive", "Unproductive",  "Unproductive"],
        'Reason': ["-", "Stopcodon", "Out-of-Frame", "Both"],
        'Percentage': [100*float(prod_reads)/total_size,
                        100*float(stopcod_inframe)/total_size,
                        100*float(outframe_nostop)/total_size,
                        100*float(both)/total_size
                        ],
        'Count': [prod_reads, stopcod_inframe, outframe_nostop, both]
    }

    _writecsv(outdir+sample_name + "_productivity2.csv", res)

def _writecsv(fname, dic):
    df = pd.DataFrame.from_dict(dic)
    df.to_csv(fname)

def process_samples(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone):
    fulldir1 = dir1+path_to_clone
    fulldir2 = dir2+path_to_clone
    fulldir3 = dir3+path_to_clone

    _process_hdf(sample_name1, fulldir1+sample_name1+"_refined_clones_annot.h5", fulldir1)
    _process_hdf(sample_name2, fulldir2+sample_name2+"_refined_clones_annot.h5", fulldir2)
    _process_hdf(sample_name3, fulldir3+sample_name3+"_refined_clones_annot.h5", fulldir3)
