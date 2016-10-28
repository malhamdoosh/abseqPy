'''
Created on 19/08/2016

@author: monther
'''
from collections import Counter

#
#  spectratype, that is, histogram of clone counts by CDR/FR nucleotide length. 
#  The spectratype is useful to detect pathological and highly clonal repertoires, 
#  as the spectratype of non-expanded T- and B-cells has a symmetric gaussian-like distribution.

def annotateSpectratypes(cloneAnnot):
    #TODO: add annotation to clonotypes, e.g., germline genes
    #TODO: add nucleotide level calculations 
    spectraTypes = {}
    # CDR1
    cdrLength = ((cloneAnnot['cdr1.end'] - cloneAnnot['cdr1.start'] + 1) / 3).astype(int)
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['cdr1'] = cdrLength
    # CDR2
    cdrLength = ((cloneAnnot['cdr2.end'] - cloneAnnot['cdr2.start'] + 1) / 3).astype(int)
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['cdr2'] = cdrLength
    #CDR3
    cdrLength = ((cloneAnnot['cdr3.end'] - cloneAnnot['cdr3.start'] + 1) / 3).astype(int)
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['cdr3'] = cdrLength
    #FR1
    cdrLength = ((cloneAnnot['fr1.end'] - cloneAnnot['fr1.start'] + 1) / 3).astype(int)
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr1'] = cdrLength
    #FR2
    cdrLength = ((cloneAnnot['fr2.end'] - cloneAnnot['fr2.start'] + 1) / 3).astype(int)
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr2'] = cdrLength
    #FR3
    cdrLength = ((cloneAnnot['fr3.end'] - cloneAnnot['fr3.start'] + 1) / 3).astype(int)
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr3'] = cdrLength
    #FR4
    cdrLength = ((cloneAnnot['fr4.end'] - cloneAnnot['fr4.start'] + 1) / 3).astype(int)
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr4'] = cdrLength
        
    return spectraTypes

# clonotype is the histogram of clone counts by CDR/FR amino acid sequence 
def annotateClonotypes(cloneSeqs):
    #TODO: add annotation to clonotypes, e.g., germline genes
    #TODO: add nucleotide level calculations
    clonoTypes = {}
    # CDR1
    seqs = cloneSeqs['cdr1'].tolist()
    clonoTypes['cdr1'] = Counter(seqs)
    # CDR2
    seqs = cloneSeqs['cdr2'].tolist()
    clonoTypes['cdr2'] = Counter(seqs)
    # CDR3
    seqs = cloneSeqs['cdr3'].tolist()
    clonoTypes['cdr3'] = Counter(seqs)
    # FR1
    seqs = cloneSeqs['fr1'].tolist()
    clonoTypes['fr1'] = Counter(seqs)
    # FR2
    seqs = cloneSeqs['fr2'].tolist()
    clonoTypes['fr2'] = Counter(seqs)
    # FR3
    seqs = cloneSeqs['fr3'].tolist()
    clonoTypes['fr3'] = Counter(seqs)
    # FR4
    seqs = cloneSeqs['fr4'].tolist()
    clonoTypes['fr4'] = Counter(seqs)
    
    
    return clonoTypes











