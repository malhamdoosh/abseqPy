'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

from abseqPy.IgRepReporting.igRepPlots import plotVenn, plotHeatmapFromDF


def generateOverlapFigures(overlapResults, noSeqs, name, siteHitsFile, stream=None):
    """

    :param overlapResults: dictionary
        dictionary with optional keys:
            {
                "order1" : {'enzyme1': {'seq_id1', 'seq_id2', 'seq_id3', ...}, 'enzyme2': {'seq_id5', ...} , ... },
                "order2" : Dataframe of n^2(all-vs-all) rows where each row is a jaccard index of the ids that each
                           pairwise comparison of the enzyme yields. This dataframe has an index column and header
                           that is identical (i.e. a "named matrix") - see calcRSAOverlapOrder2's return value
            }
            "order1" is always there, "order2" only appears if the number of enzymes is at least 3(len(sitesInfo)) >= 3)
    :param noSeqs: total number of sequences
    :param name: string. sample name
    :param siteHitsFile: string. output file name
    :param stream: logger stream
    :return: None
    """
    if overlapResults is None:
        return

    # if there's only 2 or 3 enzymes, use venn-diagram
    if "order1" in overlapResults and len(overlapResults["order1"]) in [2, 3]:
        # Ven Diagram of overlapping sequences
        title = 'Restriction sites in Sample ' + name
        title += '\nTotal is {:,}'.format(int(noSeqs))
        plotVenn(overlapResults["order1"], siteHitsFile.replace('.csv', '_venn.png'), title, stream=stream)

    # if order2 is in overlapResults, then it implies that there's AT LEAST 3 enzymes
    if "order2" in overlapResults:
        title = 'Restriction sites in Sample ' + name
        title += '\nTotal is {:,}'.format(int(noSeqs))
        # jaccard index heatmap plot
        plotHeatmapFromDF(overlapResults["order2"], siteHitsFile.replace('.csv', '_hm.png'), title=title, stream=stream)
