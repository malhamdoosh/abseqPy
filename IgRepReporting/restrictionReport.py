'''
Created on 07/02/2017

@author: monther
'''
from IgRepReporting.igRepPlots import plotVenn, plotHeatmapFromDF
from numpy import array


def generateOverlapFigures(overlapResults, noSeqs, name, siteHitsFile):
    if overlapResults is None:
        return         
    if (overlapResults.get('order1', None) is not None and 
        len(overlapResults["order1"]) in [2,3]):
        # Ven Diagram of overlapping sequences
        title = 'Restriction sites in Sample ' + name 
        title += '\nTotal is {:,}'.format(int(noSeqs)) 
        #print(stats["siteHitsSeqsIDs"])
        plotVenn(overlapResults["order1"], siteHitsFile.replace('.csv', '_venn.png'),
                 title)
    if overlapResults.get("order2", None) is not None:        
        title = 'Restriction sites in Sample ' + name 
        title += '\nTotal is {:,}'.format(int(noSeqs)) 
        #print array(overlapResults["order2"])
        plotHeatmapFromDF(overlapResults["order2"],                     
                    siteHitsFile.replace('.csv', '_hm.png'),
                    title = title) 
        
        
        
    