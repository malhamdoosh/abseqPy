'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

import sys
from IgRepertoire.igRepUtils import writeClonoTypesToFile
from IgRepReporting.igRepPlots import plotSeqLenDist, \
    generateCumulativeLogo, plotSeqDuplication, plotSeqRarefaction,\
    plotSeqRecapture, plotSeqRecaptureNew
import os
from IgRepAuxiliary.SeqUtils import createAlphabet, generateMotif
from collections import Counter

def generateDiversityReport(spectraTypes, clonoTypes, name, outDir, topClonotypes):    
    generateSpectraTypePlots(spectraTypes,  name, outDir)
    
    writeClonoTypesToFiles(clonoTypes, name, outDir, topClonotypes)
    
    estimateDiversity(clonoTypes, name, outDir)
#     generateCDRandFRLogos()
    sys.stdout.flush()

def writeClonoTypesToFiles(clonoTypes, name, outDir, topClonotypes = 100):
    print("Clonotype files are being written out ... ")
    cloneFolder = outDir + "clonotypes/"
    if (not os.path.exists(cloneFolder)):
        os.system("mkdir " + cloneFolder)
    for k in clonoTypes.keys():
        filename = cloneFolder + name + ("_%s_clonotypes_%d_over.h5" % (k, topClonotypes))
        writeClonoTypesToFile(clonoTypes[k], 
          filename, 
          topClonotypes,
          overRepresented = True)
        filename = cloneFolder + name + ("_%s_clonotypes_%d_under.h5" % (k, topClonotypes))
        writeClonoTypesToFile(clonoTypes[k], 
          filename, 
          topClonotypes,
          overRepresented = False)

def generateSpectraTypePlots(spectraTypes, name, outDir):
    specFolder = outDir + "spectratypes/"
    if (not os.path.exists(specFolder)):
        os.system("mkdir " + specFolder)
    for k in spectraTypes.keys():
        filename = specFolder + name + ('_%s_spectratype.png' % (k))
        plotSeqLenDist(spectraTypes[k], name, filename, dna=False,
              seqName=k.upper(), normed=True, maxbins=40)
        if k == 'cdr3':
            filename = specFolder + name + ('_%s_spectratype_no_outliers.png' % (k))
            plotSeqLenDist(spectraTypes[k], name, filename, dna=False,
              seqName=k.upper(), normed=True, maxbins=40, 
              removeOutliers= True)


def estimateDiversity(clonoTypes, name, outDir):
    generateSeqLogosMotifs(clonoTypes, name, outDir, "protein")
    generateRarefactionPlots(clonoTypes, name, outDir)
    calcDiversity(clonoTypes, name, outDir)

def calcDiversity(clonoTypes, name, outDir):
    print("The diversity of the library is being estimated ... ")
    regions = clonoTypes.keys()
    regions.sort()
    rarenessCounts = {}    
    for region in regions:
        rarenessCounts[region] = Counter(clonoTypes[region].values()) 
    


def generateRarefactionPlots(clonoTypes, name, outDir):
    regions = clonoTypes.keys()
    regions.sort() 
    print("Rarefaction plots are being generated .... ")
    # select CDR regions only  
    cdrWeights = []  
    cdrSeqs = []
    cdrRegions = []
    for region in regions:
        if not region.startswith("cdr"):
            continue
        cdrRegions.append(region.upper())
        cdrSeqs.append(clonoTypes[region].keys())
        cdrWeights.append(map(lambda x: clonoTypes[region][x], cdrSeqs[-1]))        
    filename = outDir + name + "_cdr_duplication.png"   
    print("\tThe duplication levels plot is being generated for CDRs .... ") 
    plotSeqDuplication(cdrWeights,                     
                     cdrRegions,
                     filename,                    
                     'Duplication of CDR Sequences')  
    print("\tThe rarefaction plot is being generated for CDRs .... ")       
    filename = outDir + name + "_cdr_rarefaction.png" 
    plotSeqRarefaction(cdrSeqs,
                     cdrRegions,
                     filename,    
                     cdrWeights,                 
                     'Rarefaction of CDR Sequences')
    print("\tThe percent recapture plot is being generated for CDRs .... ")       
    filename = outDir + name + "_cdr_recapture.png" 
    plotSeqRecaptureNew(cdrSeqs,
                     cdrRegions,
                     filename,       
                     'Percent Recapture of CDR Sequences')
    # select FR regions only
    frWeights = []  
    frSeqs = []
    frRegions = []
    for region in regions:
        if not region.startswith("fr"):
            continue
        frRegions.append(region.upper())
        frSeqs.append(clonoTypes[region].keys())
        frWeights.append(map(lambda x: clonoTypes[region][x], frSeqs[-1]))        
    filename = outDir + name + "_fr_duplication.png"   
    print("\tThe duplication levels plot is being generated for FRs .... ") 
    plotSeqDuplication(frWeights,                     
                     frRegions,
                     filename,                    
                     'Duplication of FR Sequences')  
    print("\tThe rarefaction plot is being generated for FRs .... ")       
    filename = outDir + name + "_fr_rarefaction.png" 
    plotSeqRarefaction(frSeqs,
                     frRegions,
                     filename,    
                     frWeights,                 
                     'Rarefaction of FR Sequences')
    print("\tThe percent recapture plot is being generated for FRs .... ")       
    filename = outDir + name + "_fr_recapture.png" 
    plotSeqRecaptureNew(frSeqs,
                     frRegions,
                     filename,        
                     'Percent Recapture of FR Sequences')
    # select CDR and V domain 
    cdrWeights = []  
    cdrSeqs = []
    cdrRegions = []
    for region in regions:
        if region.startswith("fr"):
            continue
        cdrRegions.append(region.upper())
        cdrSeqs.append(clonoTypes[region].keys())
        cdrWeights.append(map(lambda x: clonoTypes[region][x], cdrSeqs[-1]))        
    filename = outDir + name + "_cdr_v_duplication.png"   
    print("\tThe duplication levels plot is being generated for CDRs and V domains .... ") 
    plotSeqDuplication(cdrWeights,                     
                     cdrRegions,
                     filename,                    
                     'Duplication of CDRs and V Domains')  
    print("\tThe rarefaction plot is being generated for CDRs and V domains .... ")       
    filename = outDir + name + "_cdr_v_rarefaction.png" 
    plotSeqRarefaction(cdrSeqs,
                     cdrRegions,
                     filename,    
                     cdrWeights,                 
                     'Rarefaction of CDRs and V Domains')
    print("\tThe percent recapture plot is being generated for CDRs and V domains .... ")       
    filename = outDir + name + "_cdr_v_recapture.png" 
    plotSeqRecaptureNew(cdrSeqs,
                     cdrRegions,
                     filename,        
                     'Percent Recapture of CDRs and V Domains')
  
def generateSeqLogosMotifs(clonoTypes, name, outDir, seqType = "protein"):
    logosFolder = outDir + 'composition_logos/'    
    if (not os.path.isdir(logosFolder)):
        os.system('mkdir ' + logosFolder)
    motifsFolder =  outDir + 'motifs/'   
    if (not os.path.isdir(motifsFolder)):
        os.system('mkdir ' + motifsFolder)
    regions = clonoTypes.keys()
    regions.sort()        
    print(seqType + " sequence logos are being generated .... ")  
    for region in regions: 
        if (region == 'v'):
            continue
        print("\t" + region.upper())
        clonoType = clonoTypes[region]
        seqs = clonoType.keys()        
        weights = map(lambda x: clonoType[x], seqs)
        # Generate cumulative sequence logos using Toby's approach
        #TODO: generate composition logos by IGV family
        filename = logosFolder + name + ("_%s_cumulative_logo.png" % (region))        
        generateCumulativeLogo(seqs, weights, region, filename)
        # Generate sequence motif logos using weblogo                
        # generate logos without alignment
        filename = motifsFolder + name + ("_%s_motif_logo.png" % (region))
        alphabet = createAlphabet(align=False, protein=True, extendAlphabet = True)
        m = generateMotif(seqs, region, alphabet, filename,  align = False,
                          protein = True, weights= weights)
        # generate  logos after alignment
        filename = motifsFolder + name + ("_%s_motif_aligned_logo.png" % (region))
        alphabet = createAlphabet(align=True, protein=True, extendAlphabet = True)
        m = generateMotif(seqs, region, alphabet, filename,  align = True,
                          protein = True, weights= weights)
    
    

  
# quantify CDR sequence diversity      
          
#         if (not exists(self.outputDir + self.name + 
#                  '_Vdomain_diversity.png')):            
#     #         i = 0
#             VH = []
#             for (id, f1, c1, f2, c2, f3, c3, f4) in zip(self.cloneSeqs.index.tolist(),
#                                                         self.cloneSeqs['fr1'].tolist(),
#                                                           self.cloneSeqs['cdr1'].tolist(),
#                                                           self.cloneSeqs['fr2'].tolist(),
#                                                           self.cloneSeqs['cdr2'].tolist(),
#                                                           self.cloneSeqs['fr3'].tolist(),
#                                                           self.cloneSeqs['cdr3'].tolist(),
#                                                           self.cloneSeqs['fr4'].tolist()):           
#                 try:
#                     VH += [''.join([f1, c1, f2, c2, f3, c3, f4])]
#                 except:
#                     if (f4 is None or isnan(f4)):  # or c3 is None or isnan(c3):
#                         VH += [''.join([f1, c1, f2, c2, f3, c3])]
#                     else:
#                         print(id, f1, c1, f2, c2, f3, c3, f4)
# #                 i += 1
# #         print(i)
# #         sys.exit()
#             plotSeqDuplication([self.cloneSeqs['cdr1'].tolist(),
#                               self.cloneSeqs['cdr2'].tolist(),
#                               self.cloneSeqs['cdr3'].tolist(),
#                               VH],
#                              self.outputDir + self.name + 
#                      '_Vdomain__Vdomain_ication.png',
#                              ['CDR1', 'CDR2', 'CDR3', 'V Domain'],
#                              'Duplication of V Domain Sequences')
#             plotSeqRarefaction([self.cloneSeqs['cdr1'].tolist(),
#                           self.cloneSeqs['cdr2'].tolist(),
#                           self.cloneSeqs['cdr3'].tolist(),
#                           VH],
#                          self.outputDir + self.name + 
#                  '_Vdomain_diversity.png',
#                          ['CDR1', 'CDR2', 'CDR3', 'V Domain'],
#                          'Diversity of V Domain Sequences')
#         gc.collect()
#         
#         plotSeqDuplication([self.cloneSeqs['cdr1'].tolist(),
#                           self.cloneSeqs['cdr2'].tolist(),
#                           self.cloneSeqs['cdr3'].tolist()
#                           ],
#                          self.outputDir + self.name + 
#                  '_cdr_duplication.png',
#                          ['CDR1', 'CDR2', 'CDR3'],
#                          'Duplication of CDR Sequences')        
#         
#         plotSeqRarefaction([self.cloneSeqs['cdr1'].tolist(),
#                           self.cloneSeqs['cdr2'].tolist(),
#                           self.cloneSeqs['cdr3'].tolist()
#                         ],
#                          self.outputDir + self.name + 
#                  '_cdr_diversity.png',
#                          ['CDR1', 'CDR2', 'CDR3'],
#                          'Diversity of CDR Sequences')
#         gc.collect()
   
   
# Quantify FR sequence diversity
#         plotSeqDuplication([self.cloneSeqs['fr1'].tolist(),
#                           self.cloneSeqs['fr2'].tolist(),
#                           self.cloneSeqs['fr3'].tolist(),
#                           self.cloneSeqs['fr4'].tolist()],
#                          self.outputDir + self.name + 
#                  '_fr_duplication.png',
#                          ['FR1', 'FR2', 'FR3', 'FR4'],
#                          'Duplication of FR Sequences')
#         gc.collect()
#         plotSeqRarefaction([self.cloneSeqs['fr1'].tolist(),
#                           self.cloneSeqs['fr2'].tolist(),
#                           self.cloneSeqs['fr3'].tolist(),
#                           self.cloneSeqs['fr4'].tolist()],
#                          self.outputDir + self.name + 
#                  '_fr_diversity.png',
#                          ['FR1', 'FR2', 'FR3', 'FR4'],
#                          'Diversity of FR Sequences')


# quantify V domain sequence diversity                
#         if (not exists(self.outputDir + self.name + 
#                  '_Vdomain_duplication_family.png')):
#             print("Grouping V domain sequences per family ...")
#             VH = {}
#     #         i = 0
#             ighvs = map(lambda x : x.split('-')[0].split('/')[0], self.cloneSeqs['germline'].tolist())
#             for ighv in set(ighvs):
#                 VH[ighv] = []
#             for (ighv, f1, c1, f2, c2, f3, c3, f4) in zip(ighvs,
#                                                         self.cloneSeqs['fr1'].tolist(),
#                                                           self.cloneSeqs['cdr1'].tolist(),
#                                                           self.cloneSeqs['fr2'].tolist(),
#                                                           self.cloneSeqs['cdr2'].tolist(),
#                                                           self.cloneSeqs['fr3'].tolist(),
#                                                           self.cloneSeqs['cdr3'].tolist(),
#                                                           self.cloneSeqs['fr4'].tolist()):           
#                 try:
#                     VH[ighv].append(''.join([f1, c1, f2, c2, f3, c3, f4]))
#                 except:
#                     if (f4 is None or isnan(f4)):  # or c3 is None or isnan(c3):
#                         VH += [''.join([f1, c1, f2, c2, f3, c3])]
#                     else:
#                         print(id, f1, c1, f2, c2, f3, c3, f4)
#             ighvs = VH.keys()
#             ighvs.sort()        
#             
#             plotSeqDuplication(map(lambda x:VH[x], ighvs),
#                              self.outputDir + self.name + 
#                      '_Vdomain_duplication_family.png',
#                              ighvs,
#                              'Duplication of V Domain Sequences Per Family', True)
#             gc.collect()