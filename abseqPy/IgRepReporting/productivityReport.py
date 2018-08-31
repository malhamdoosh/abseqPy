'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

import gc
import pandas as pd
import os

from collections import Counter, OrderedDict
from numpy import nan

from abseqPy.IgRepReporting.igRepPlots import plotDist
from abseqPy.IgRepertoire.igRepUtils import compressCountsFamilyLevel


def generateProductivityReport(cloneAnnot, cloneSeqs, name, chain, outputDir, stream=None):
    # since np.nan is considered different objects, canonicalize them using 'NaN' string representation
    nanString = 'NaN'
    cloneAnnot.fillna(nanString, inplace=True)

    productive = extractProductiveClones(cloneAnnot, name, outputDir, stream=stream)
    productiveFamilyDist = compressCountsFamilyLevel(Counter(productive['vgene'].tolist()))
    plotDist(productiveFamilyDist, name, os.path.join(outputDir, name + '_igv_dist_productive.csv'),
             title='IGV Abundance of Productive Clones',
             proportion=True, stream=stream)
    del productiveFamilyDist
    writeProdStats(cloneAnnot, name, outputDir)
    writeCDRStats(productive, name, outputDir, suffix='productive', stream=stream)
    writeFRStats(productive, name, outputDir, suffix='productive', stream=stream)
    writeGeneStats(productive, name, chain, outputDir, suffix='productive', stream=stream)
    writeStopCodonStats(cloneAnnot, cloneSeqs, name, outputDir, inframe=True, stream=stream)
    writeStopCodonStats(cloneAnnot, cloneSeqs, name, outputDir, inframe=False, stream=stream)

    # now that counting is complete, replace all 'NaN' strings with np.nan again
    cloneAnnot.replace(nanString, nan, inplace=True)


def writeProdStats(cloneAnnot, sampleName, outdir):
    """
    writes the statistics of un-productiveness vs productiveness.Gives the number of productive clones
    vs unproductive clones with "reason" for being unproductive - i.e. out-of-frame/stopcodon/both
    :param cloneAnnot: refined_clones_annot.h5's dataframe object
    :param sampleName: name of this sample being analysed
    :param outdir: output directory
    :return: None. Produces a csv file in outdir
    """
    stopcod_inframe = len(cloneAnnot[(cloneAnnot['v-jframe'] == 'In-frame') & (cloneAnnot['stopcodon'] == 'Yes')])
    outframe_nostop = len(cloneAnnot[(cloneAnnot['v-jframe'] == 'Out-of-frame') & (cloneAnnot['stopcodon'] == 'No')])
    both = len(cloneAnnot[(cloneAnnot['v-jframe'] == 'Out-of-frame') & (cloneAnnot['stopcodon'] == 'Yes')])
    prod_reads = len(cloneAnnot[(cloneAnnot['v-jframe'] == 'In-frame') & (cloneAnnot['stopcodon'] == 'No')])

    # percentage calculated from total sample size
    total_size = len(cloneAnnot)
    if total_size:
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

        df = pd.DataFrame.from_dict(res)
        df.to_csv(os.path.join(outdir, sampleName + "_productivity.csv"))


def writeGeneStats(cloneAnnot, name, chain, outputDir, suffix, stream=None):
    # V gene stats
    gaps = Counter(cloneAnnot['vgaps'].tolist())
    plotDist(gaps, name, os.path.join(outputDir, name +
             '_igv_gaps_dist.csv'), title='Gaps in V Gene',
             proportion=True, rotateLabels=False, top=20, stream=stream)
    mismatches = Counter(cloneAnnot['vmismatches'].tolist())
    plotDist(mismatches, name, os.path.join(outputDir, name +
             '_igv_mismatches_dist.csv'), title='Mismatches in V Gene',
             proportion=True, rotateLabels=False, top=20, stream=stream)
    # D gene stats
    if chain == 'hv':
        gaps = Counter(cloneAnnot['dgaps'].tolist())
        plotDist(gaps, name, os.path.join(outputDir, name +
                 '_igd_gaps_dist.csv'), title='Gaps in D Gene',
                 proportion=False, rotateLabels=False, stream=stream)
        mismatches = Counter(cloneAnnot['dmismatches'].tolist())
#         print(mismatches)
        plotDist(mismatches, name, os.path.join(outputDir, name +
                 '_igd_mismatches_dist.csv'), title='Mismatches in D Gene',
                 proportion=False, rotateLabels=False, stream=stream)
    # J gene stats
    gaps = Counter(cloneAnnot['jgaps'].tolist())
    plotDist(gaps, name, os.path.join(outputDir, name +
             '_igj_gaps_dist.csv'), title='Gaps in J Gene',
             proportion=False, rotateLabels=False, stream=stream)
    mismatches = Counter(cloneAnnot['jmismatches'].tolist())
    plotDist(mismatches, name, os.path.join(outputDir, name +
             '_igj_mismatches_dist.csv'), title='Mismatches in J Gene',
             proportion=False, rotateLabels=False, stream=stream)


def writeCDRStats(cloneAnnot, name, outputDir, suffix = '', stream=None):
    # CDR1 statistics
    cdrGaps = Counter(cloneAnnot['cdr1.gaps'].tolist())
    plotDist(cdrGaps, name, os.path.join(outputDir, name +
             '_cdr1_gaps_dist.csv'), title='Gaps in CDR1',
             proportion=False, rotateLabels=False, stream=stream)
    cdrMismatches = Counter(cloneAnnot['cdr1.mismatches'].tolist())
    plotDist(cdrMismatches, name, os.path.join(outputDir, name +
             '_cdr1_mismatches_dist.csv'), title='Mismatches in CDR1',
             proportion=False, rotateLabels=False, stream=stream)
    # CDR2 stats    
    cdrGaps = Counter(cloneAnnot['cdr2.gaps'].tolist())
    plotDist(cdrGaps, name, os.path.join(outputDir, name +
             '_cdr2_gaps_dist.csv'), title='Gaps in CDR2',
             proportion=False, rotateLabels=False, stream=stream)
    cdrMismatches = Counter(cloneAnnot['cdr2.mismatches'].tolist())
    plotDist(cdrMismatches, name, os.path.join(outputDir, name +
             '_cdr2_mismatches_dist.csv'), title='Mismatches in CDR2',
             proportion=False, rotateLabels=False, stream=stream)
    # CDR3 stats
    cdrGaps = Counter(cloneAnnot['cdr3g.gaps'])
#         print(len(cdrGaps))
    plotDist(cdrGaps, name, os.path.join(outputDir, name +
             '_cdr3_gaps_dist.csv'), title='Gaps in CDR3 (Germline)',
             proportion=False, rotateLabels=False, stream=stream)
    cdrMismatches = Counter(cloneAnnot['cdr3g.mismatches'].tolist())
    plotDist(cdrMismatches, name, os.path.join(outputDir, name +
             '_cdr3_mismatches_dist.csv'), title='Mismatches in CDR3 (Germline)',
             proportion=False, rotateLabels=False, stream=stream)
    gc.collect()


def writeFRStats(cloneAnnot, name, outputDir, suffix = '', stream=None):
    # FR1 statistics 
    gaps = Counter(cloneAnnot['fr1.gaps'].tolist())
    plotDist(gaps, name, os.path.join(outputDir,  name +
             '_fr1_gaps_dist.csv'), title='Gaps in FR1',
             proportion=False, rotateLabels=False, stream=stream)
    mismatches = Counter(cloneAnnot['fr1.mismatches'].tolist())
    plotDist(mismatches, name, os.path.join(outputDir, name +
             '_fr1_mismatches_dist.csv'), title='Mismatches in FR1',
             proportion=False, rotateLabels=False, stream=stream)
    # FR2 statistics 
    gaps = Counter(cloneAnnot['fr2.gaps'].tolist())
    plotDist(gaps, name, os.path.join(outputDir, name +
             '_fr2_gaps_dist.csv'), title='Gaps in FR2',
             proportion=False, rotateLabels=False, stream=stream)
    mismatches = Counter(cloneAnnot['fr2.mismatches'].tolist())
    plotDist(mismatches, name, os.path.join(outputDir, name +
             '_fr2_mismatches_dist.csv'), title='Mismatches in FR2',
             proportion=False, rotateLabels=False, stream=stream)
    # FR3 statistics 
    gaps = Counter(cloneAnnot['fr3g.gaps'].tolist())
    plotDist(gaps, name, os.path.join(outputDir, name +
             '_fr3_gaps_dist.csv'), title='Gaps in FR3 (Germline)',
             proportion=False, rotateLabels=False, stream=stream)
    mismatches = Counter(cloneAnnot['fr3g.mismatches'].tolist())
    plotDist(mismatches, name, os.path.join(outputDir, name +
             '_fr3_mismatches_dist.csv'), title='Mismatches in FR3 (Germline)',
             proportion=False, rotateLabels=False, stream=stream)
    
    gc.collect()   


def extractProductiveClones(cloneAnnot, name, outputDir, stream=None):
    # v-j rearrangement frame distribution 
    vjframeDist = Counter(cloneAnnot['v-jframe'].tolist())        
    plotDist(vjframeDist, name, os.path.join(outputDir, name +
             '_vjframe_dist.csv'), title='V-D-J Rearrangement',
             proportion=False, rotateLabels=False, stream=stream)
    del vjframeDist
    # plot the family distribution of out-of-frame
    outOfFrame = cloneAnnot[cloneAnnot['v-jframe'] != 'In-frame']
    outOfFrameFamilyDist = compressCountsFamilyLevel(Counter(outOfFrame['vgene'].tolist()))
    plotDist(outOfFrameFamilyDist, name, os.path.join(outputDir, name +
             '_igv_dist_out_of_frame.csv'),
              title='IGV Abundance of Out-Of-frame Clones',
             proportion=True, stream=stream)
    del outOfFrameFamilyDist
    # Indels in CDR1 and FR1    
    cdrGaps = Counter(outOfFrame['cdr1.gaps'].tolist())
    plotDist(cdrGaps, name, os.path.join(outputDir, name +
             '_cdr1_gaps_dist_out_of_frame.csv'), title='Gaps in CDR1',
             proportion=False, rotateLabels=False, stream=stream)
    frGaps = Counter(outOfFrame['fr1.gaps'].tolist())
    plotDist(frGaps, name, os.path.join(outputDir, name +
             '_fr1_gaps_dist_out_of_frame.csv'), title='Gaps in FR1',
             proportion=False, rotateLabels=False, stream=stream)
    del  cdrGaps, frGaps
    # Indels in CDR2 and FR2
    cdrGaps = Counter(outOfFrame['cdr2.gaps'].tolist())
    plotDist(cdrGaps, name, os.path.join(outputDir, name +
             '_cdr2_gaps_dist_out_of_frame.csv'), title='Gaps in CDR2',
             proportion=False, rotateLabels=False, stream=stream)
    frGaps = Counter(outOfFrame['fr2.gaps'].tolist())
    plotDist(frGaps, name, os.path.join(outputDir, name +
             '_fr2_gaps_dist_out_of_frame.csv'), title='Gaps in FR2',
             proportion=False, rotateLabels=False, stream=stream)
    del cdrGaps, frGaps
    # Indels in CDR3 and FR3
    cdrGaps = Counter(outOfFrame['cdr3g.gaps'])
#         print(len(cdrGaps))
    plotDist(cdrGaps, name, os.path.join(outputDir, name +
             '_cdr3_gaps_dist_out_of_frame.csv'), title='Gaps in CDR3 (Germline)',
             proportion=False, rotateLabels=False, stream=stream)
    frGaps = Counter(outOfFrame['fr3g.gaps'].tolist())
    plotDist(frGaps, name, os.path.join(outputDir, name +
             '_fr3_gaps_dist_out_of_frame.csv'), title='Gaps in FR3 (Germline)',
             proportion=False, rotateLabels=False, stream=stream)
    del cdrGaps, frGaps
#     # Indels in FR4
#     frGaps = Counter(outOfFrame['fr3.gaps'].tolist())
#     plotDist(frGaps, name, outputDir + name + 
#              '_fr3_gaps_dist_out_of_frame.csv', title='Gaps in FR3',
#              proportion=False, rotateLabels=False)   
    del outOfFrame    
    # choose only In-frame RNA clones
    inFrame = cloneAnnot[cloneAnnot['v-jframe'] == 'In-frame']
    # Stop Codon 
    stopcodonInFrameDist = Counter(inFrame['stopcodon'].tolist())
    plotDist(stopcodonInFrameDist, name, os.path.join(outputDir,  name +
             '_stopcodon_dist_in_frame.csv'), title='Stop Codons in In-frame Clones',
             proportion=False, rotateLabels=False, stream=stream)
    
    # stop codon family distribution
    stopcodFamily = Counter(inFrame[inFrame['stopcodon'] == 'Yes']['vgene'].tolist())
    stopcodFamily = compressCountsFamilyLevel(stopcodFamily)
    plotDist(stopcodFamily, name, os.path.join(outputDir, name +
             '_igv_dist_inframe_unproductive.csv'),
             title='IGV Abundance of In-frame Unproductive Clones',
             proportion=True, stream=stream)
    del stopcodonInFrameDist, stopcodFamily
#         print(stopcodFamily)
    # choose only productive RNA sequences 
    productive = inFrame[inFrame['stopcodon'] == 'No']
    gc.collect()
    
    return productive
    
  
def writeStopCodonStats(cloneAnnot, cloneSeqs, name, outputDir, inframe, stream=None):
    """
    This function maintains the hypothesis that a stop codon is independent of
    previous stop codons. It increments the counter for each region as long as there's
    AT LEAST ONE stop codon in the specified region. This is especially true if the sequence
    is in-frame.
    :param cloneAnnot: .*_clone_annot.h5
    :param cloneSeqs: .*_clones_seq.h5
    :param name: sample name
    :param outputDir: output directory
    :param inframe: True if only for inframe sequences, false if only for out-of-frame sequences
    :param stream: debugging stream
    :return:
    """
    regions = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']

    counter = {}
    frameStatus = 'In-frame' if inframe else 'Out-of-frame'
    cloneSeqs = cloneSeqs.loc[cloneAnnot[cloneAnnot['v-jframe'] == frameStatus].index]
    for region in regions:
        counter[region] = sum(cloneSeqs[region.lower()].str.contains("*", regex=False))
    orderedCounter = OrderedDict((reg, counter[reg]) for reg in regions)
    plotDist(orderedCounter, name, os.path.join(outputDir, name
             + '_stopcodon_region_{}.csv').format('inframe' if inframe else 'outframe'),
             title="Stop codon in FRs and CDRs of {} sequences".format(frameStatus),
             proportion=True, sortValues=False, maintainx=True, stream=stream)
