'''
Created on 15/08/2016

@author: monther
'''
import sys


def generateDiversityReport(cloneAnnot, cloneSeqs):
    generateSpectraTypes()
    
#     generateCDRandFRLogos()
    sys.stdout.flush()


#
#  spectratype, that is, histogram of read counts by CDR3 nucleotide length. 
#  The spectratype is useful to detect pathological and highly clonal repertoires, 
#  as the spectratype of non-expanded T- and B-cells has a symmetric gaussian-like distribution.

def generateSpectraTypes(cloneAnnot, cloneSeqs):
    # CDR1
    cdrLength = (self.cloneAnnot['cdr1.end'] - self.cloneAnnot['cdr1.start'] + 1) / 3
    cdrLength = cdrLength.tolist()
    histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
             '_cdr1_len_dist.png', dna=False,
              seqName='CDR1', normed=True, maxbins=20)
    # CDR2 stats
    cdrLength = (self.cloneAnnot['cdr2.end'] - self.cloneAnnot['cdr2.start'] + 1) / 3
    cdrLength = cdrLength.tolist()
    histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
             '_cdr2_len_dist.png', dna=False,
             seqName='CDR2', normed=True, maxbins=20)
    # CDR3 stats
    cdrLength = (self.cloneAnnot['cdr3.end'] - self.cloneAnnot['cdr3.start'] + 1) / 3
    cdrLength = cdrLength.tolist()
    histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
             '_cdr3_len_dist_nooutliers.png', dna=False,
             seqName='CDR3', normed=True, maxbins=20,
             removeOutliers=True)
    histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name + 
             '_cdr3_len_dist.png', dna=False,
             seqName='CDR3', normed=True, maxbins=20,
             removeOutliers=False)
    gc.collect()
    sys.stdout.flush()
    # FR1 statistics 
    frLength = (self.cloneAnnot['fr1.end'] - self.cloneAnnot['fr1.start'] + 1) / 3
    frLength = frLength.tolist()
    histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
             '_fr1_len_dist.png', dna=False,
              seqName='FR1', normed=True, maxbins=20)
    # FR2 statistics 
    frLength = (self.cloneAnnot['fr2.end'] - self.cloneAnnot['fr2.start'] + 1) / 3
    frLength = frLength.tolist()
    histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
             '_fr2_len_dist.png', dna=False,
              seqName='FR2', normed=True, maxbins=20)
    # FR3 statistics 
    frLength = (self.cloneAnnot['fr3.end'] - self.cloneAnnot['fr3.start'] + 1) / 3
    frLength = frLength.tolist()
    histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
             '_fr3_len_dist.png', dna=False,
              seqName='FR3', normed=True, maxbins=20)
    # FR4
    frLength = (self.cloneAnnot['fr4.end'] - self.cloneAnnot['fr4.start'] + 1) / 3
    frLength = frLength.tolist()
    histcals = plotSeqLenDist(frLength, self.name, self.outputDir + self.name + 
             '_fr4_len_dist.png', dna=False,
              seqName='FR4', normed=True, maxbins=20)
    sys.stdout.flush()
    
# def generateCDRandFRLogos(self):
# #         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'  
# #         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
#     seqs = {}
#     for k in self.cloneSeqs.columns:
#         if k.startswith('cdr') :
#             seqs[k] = self.cloneSeqs[k].tolist()
#         if k.startswith('fr'):
#             seqs[k] = self.cloneSeqs[k].tolist()
#     # generate Toby's logos ?!?!
#     generateProteinLogos(seqs, self.outputDir + self.name + '_cdr_fr_Toby')
#     
#     # generate CDR/FR logos without alignment
#     generateMotifs(seqs, False, self.outputDir + self.name + '_cdr_fr',
#                    protein=True) 
#     # generate CDR/FR logos after alignment
#     generateMotifs(seqs, True,
#                 self.outputDir + self.name + '_cdr_fr_aligned',
#                 protein=True)
#     # generate CDR3 logos per  germline 
#     seqs = {}
#     vgenes = self.cloneSeqs["germline"].tolist()
#     vgenes = map(lambda x: x.split('*')[0], vgenes)
#     for vgene in set(vgenes):
#         seqs["cdr3_"+vgene.replace("/", "-")] = self.cloneSeqs.loc[map(lambda x: x == vgene, vgenes), "cdr3"].tolist()            
#     # generate Toby's logos ?!?!
#     generateProteinLogos(seqs, self.outputDir + self.name + '_cdr_Toby_gene')
#     # generate CDR/FR logos after alignment
#     generateMotifs(seqs, True,
#                 self.outputDir + self.name + '_cdr_gene_aligned',
#                 protein=True)
#     # generate CDR3 logos per family        
#     seqs = {}
#     vfams = map(lambda x: x.split('-')[0].split('/')[0], vgenes)
#     for vfam in set(vfams):
#         seqs["cdr3_"+vfam] = self.cloneSeqs.loc[map(lambda x: x == vgene, vgenes), "cdr3"].tolist()
#     # generate Toby's logos ?!?!
#     generateProteinLogos(seqs, self.outputDir + self.name + '_cdr_fr_Toby')
#     # generate CDR/FR logos after alignment
#     generateMotifs(seqs, True,
#                 self.outputDir + self.name + '_cdr_fr_aligned',
#                 protein=True)
#        
    

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
#             plotSeqDiversity([self.cloneSeqs['cdr1'].tolist(),
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
#         plotSeqDiversity([self.cloneSeqs['cdr1'].tolist(),
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
#         plotSeqDiversity([self.cloneSeqs['fr1'].tolist(),
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