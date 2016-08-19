'''
Created on 05/08/2016

@author: monther
'''
from collections import Counter
from IgRepReporting.igRepPlots import plotDist
from IgRepertoire.igRepUtils import compressCountsFamilyLevel
from numpy import isnan
import gc


def generateProductivityReport(cloneAnnot, name, outputDir):
    productive = extractProductiveClones(cloneAnnot, name, outputDir)
    productiveFamilyDist = compressCountsFamilyLevel(Counter(productive['vgene'].tolist()))
    plotDist(productiveFamilyDist, name, outputDir + name + 
             '_igv_dist_productive.png',
              title='IGV Abundance of Productive Clones',
             proportion=True)
    del productiveFamilyDist
    writeCDRStats(productive, name, outputDir, suffix = 'productive')
    writeFRStats(productive, name, outputDir, suffix = 'productive')
    writeGeneStats(productive, name, outputDir, suffix = 'productive')
    
def writeGeneStats(cloneAnnot, name, outputDir, suffix):
    # V gene stats
    gaps = Counter(cloneAnnot['vgaps'].tolist())
    plotDist(gaps, name, outputDir + name + 
             '_igv_gaps_dist.png', title='Gaps in V Gene',
             proportion=True, rotateLabels=False, top=20) 
    mismatches = Counter(cloneAnnot['vmismatches'].tolist())
    plotDist(mismatches, name, outputDir + name + 
             '_igv_mismatches_dist.png', title='Mismatches in V Gene',
             proportion=True, rotateLabels=False, top=20) 
    # D gene stats
    gaps = Counter(cloneAnnot['dgaps'].tolist())
    plotDist(gaps, name, outputDir + name + 
             '_igd_gaps_dist.png', title='Gaps in D Gene',
             proportion=False, rotateLabels=False) 
    mismatches = Counter(cloneAnnot['dmismatches'].tolist())
    plotDist(mismatches, name, outputDir + name + 
             '_igd_mismatches_dist.png', title='Mismatches in D Gene',
             proportion=False, rotateLabels=False)
    # J gene stats
    gaps = Counter(cloneAnnot['jgaps'].tolist())
    plotDist(gaps, name, outputDir + name + 
             '_igj_gaps_dist.png', title='Gaps in J Gene',
             proportion=False, rotateLabels=False) 
    mismatches = Counter(cloneAnnot['jmismatches'].tolist())
    plotDist(mismatches, name, outputDir + name + 
             '_igj_mismatches_dist.png', title='Mismatches in J Gene',
             proportion=False, rotateLabels=False)
  
def writeCDRStats(cloneAnnot, name, outputDir, suffix = ''):
    # CDR1 statistics
    cdrGaps = Counter(cloneAnnot['cdr1.gaps'].tolist())
    plotDist(cdrGaps, name, outputDir + name + 
             '_cdr1_gaps_dist.png', title='Gaps in CDR1',
             proportion=False, rotateLabels=False) 
    cdrMismatches = Counter(cloneAnnot['cdr1.mismatches'].tolist())
    plotDist(cdrMismatches, name, outputDir + name + 
             '_cdr1_mismatches_dist.png', title='Mismatches in CDR1',
             proportion=False, rotateLabels=False) 
    # CDR2 stats    
    cdrGaps = Counter(cloneAnnot['cdr2.gaps'].tolist())
    plotDist(cdrGaps, name, outputDir + name + 
             '_cdr2_gaps_dist.png', title='Gaps in CDR2',
             proportion=False, rotateLabels=False)
    cdrMismatches = Counter(cloneAnnot['cdr2.mismatches'].tolist())
    plotDist(cdrMismatches, name, outputDir + name + 
             '_cdr2_mismatches_dist.png', title='Mismatches in CDR2',
             proportion=False, rotateLabels=False)
    # CDR3 stats
    cdrGaps = Counter([x if not isnan(x) else 'NA' for x in cloneAnnot['cdr3.gaps'] ])
#         print(len(cdrGaps))
    plotDist(cdrGaps, name, outputDir + name + 
             '_cdr3_gaps_dist.png', title='Gaps in CDR3 (Germline)',
             proportion=False, rotateLabels=False)
    cdrMismatches = Counter(cloneAnnot['cdr3.mismatches'].tolist())
    plotDist(cdrMismatches, name, outputDir + name + 
             '_cdr3_mismatches_dist.png', title='Mismatches in CDR3 (Germline)',
             proportion=False, rotateLabels=False)
    gc.collect()
 
def writeFRStats(cloneAnnot, name, outputDir, suffix = ''):
    # FR1 statistics 
    gaps = Counter(cloneAnnot['fr1.gaps'].tolist())
    plotDist(gaps, name, outputDir + name + 
             '_fr1_gaps_dist.png', title='Gaps in FR1',
             proportion=False, rotateLabels=False) 
    mismatches = Counter(cloneAnnot['fr1.mismatches'].tolist())
    plotDist(mismatches, name, outputDir + name + 
             '_fr1_mismatches_dist.png', title='Mismatches in FR1',
             proportion=False, rotateLabels=False) 
    # FR2 statistics 
    gaps = Counter(cloneAnnot['fr2.gaps'].tolist())
    plotDist(gaps, name, outputDir + name + 
             '_fr2_gaps_dist.png', title='Gaps in FR2',
             proportion=False, rotateLabels=False) 
    mismatches = Counter(cloneAnnot['fr2.mismatches'].tolist())
    plotDist(mismatches, name, outputDir + name + 
             '_fr2_mismatches_dist.png', title='Mismatches in FR2',
             proportion=False, rotateLabels=False) 
    # FR3 statistics 
    gaps = Counter(cloneAnnot['fr3.gaps'].tolist())
    plotDist(gaps, name, outputDir + name + 
             '_fr3_gaps_dist.png', title='Gaps in FR3',
             proportion=False, rotateLabels=False) 
    mismatches = Counter(cloneAnnot['fr3.mismatches'].tolist())
    plotDist(mismatches, name, outputDir + name + 
             '_fr3_mismatches_dist.png', title='Mismatches in FR3',
             proportion=False, rotateLabels=False)
    
    gc.collect()   

def extractProductiveClones(cloneAnnot, name, outputDir):
    # v-j rearrangement frame distribution 
    vjframeDist = Counter(cloneAnnot['v-jframe'].tolist())        
    plotDist(vjframeDist, name, outputDir + name + 
             '_vjframe_dist.png', title='V-D-J Rearrangement',
             proportion=False, rotateLabels=False)    
    del vjframeDist
    # plot the family distribution of out-of-frame
    outOfFrame = cloneAnnot[cloneAnnot['v-jframe'] != 'In-frame']
    outOfFrameFamilyDist = compressCountsFamilyLevel(Counter(outOfFrame['vgene'].tolist()))
    plotDist(outOfFrameFamilyDist, name, outputDir + name + 
             '_igv_dist_out_of_frame.png',
              title='IGV Abundance of Out-Of-frame Clones',
             proportion=True)
    del outOfFrameFamilyDist
    # Indels in CDR1 and FR1    
    cdrGaps = Counter(outOfFrame['cdr1.gaps'].tolist())
    plotDist(cdrGaps, name, outputDir + name + 
             '_cdr1_gaps_dist_out_of_frame.png', title='Gaps in CDR1',
             proportion=False, rotateLabels=False)
    frGaps = Counter(outOfFrame['fr1.gaps'].tolist())
    plotDist(frGaps, name, outputDir + name + 
             '_fr1_gaps_dist_out_of_frame.png', title='Gaps in FR1',
             proportion=False, rotateLabels=False)
    del  cdrGaps, frGaps
    # Indels in CDR2 and FR2
    cdrGaps = Counter(outOfFrame['cdr2.gaps'].tolist())
    plotDist(cdrGaps, name, outputDir + name + 
             '_cdr2_gaps_dist_out_of_frame.png', title='Gaps in CDR2',
             proportion=False, rotateLabels=False)
    frGaps = Counter(outOfFrame['fr2.gaps'].tolist())
    plotDist(frGaps, name, outputDir + name + 
             '_outframe_fr2_gaps_dist.png', title='Gaps in FR2',
             proportion=False, rotateLabels=False)
    del cdrGaps, frGaps
    # Indels in CDR3 and FR3
    cdrGaps = Counter([x if not isnan(x) else 'NA' for x in outOfFrame['cdr3.gaps'] ])
#         print(len(cdrGaps))
    plotDist(cdrGaps, name, outputDir + name + 
             '_cdr3_gaps_dist_out_of_frame.png', title='Gaps in CDR3 (Germline)',
             proportion=False, rotateLabels=False)
    frGaps = Counter(outOfFrame['fr3.gaps'].tolist())
    plotDist(frGaps, name, outputDir + name + 
             '_fr3_gaps_dist_out_of_frame.png', title='Gaps in FR3',
             proportion=False, rotateLabels=False)
    del cdrGaps, frGaps
#     # Indels in FR4
#     frGaps = Counter(outOfFrame['fr3.gaps'].tolist())
#     plotDist(frGaps, name, outputDir + name + 
#              '_fr3_gaps_dist_out_of_frame.png', title='Gaps in FR3',
#              proportion=False, rotateLabels=False)   
    del outOfFrame    
    # choose only In-frame RNA clones
    inFrame = cloneAnnot[cloneAnnot['v-jframe'] == 'In-frame']
    # Stop Codon 
    stopcodonInFrameDist = Counter(inFrame['stopcodon'].tolist())
    plotDist(stopcodonInFrameDist, name, outputDir + name + 
             '_stopcodon_dist_in_frame.png', title='Stop Codons in In-frame Clones',
             proportion=False, rotateLabels=False)
    
    # stop codon family distribution
    stopcodFamily = Counter(inFrame[inFrame['stopcodon'] == 'Yes']['vgene'].tolist())
    stopcodFamily = compressCountsFamilyLevel(stopcodFamily)
    plotDist(stopcodFamily, name, outputDir + name + 
             '_igv_dist_inframe_unproductive.png',
              title='IGV Abundance of In-frame Unproductive Clones',
             proportion=True)
    del stopcodonInFrameDist, stopcodFamily
#         print(stopcodFamily)
    # choose only productive RNA sequences 
    productive = inFrame[inFrame['stopcodon'] == 'No']
    gc.collect()
    
    return productive
    
  


