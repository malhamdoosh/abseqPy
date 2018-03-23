'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

from numpy import array

from abseq.IgRepReporting.igRepPlots import plotVenn, plotHeatmapFromDF


def generateOverlapFigures(overlapResults, noSeqs, name, siteHitsFile, stream=None):
    if overlapResults is None:
        return
    if overlapResults.get('order1', None) is not None and len(overlapResults["order1"]) in [2, 3]:
        # Ven Diagram of overlapping sequences
        title = 'Restriction sites in Sample ' + name
        title += '\nTotal is {:,}'.format(int(noSeqs))
        # print(stats["siteHitsSeqsIDs"])
        plotVenn(overlapResults["order1"], siteHitsFile.replace('.csv', '_venn.png'), title, stream=stream)
    if overlapResults.get("order2", None) is not None:
        title = 'Restriction sites in Sample ' + name
        title += '\nTotal is {:,}'.format(int(noSeqs))
        # print array(overlapResults["order2"])
        plotHeatmapFromDF(overlapResults["order2"], siteHitsFile.replace('.csv', '_hm.png'), title=title, stream=stream)
