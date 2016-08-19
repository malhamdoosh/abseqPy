'''
Created on 19/08/2016

@author: monther
'''
from collections import Counter

def identifySpectratypes(cloneAnnot):
    spectraTypes = {}
    # CDR1
    cdrLength = (cloneAnnot['cdr1.end'] - cloneAnnot['cdr1.start'] + 1) / 3
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['cdr1'] = cdrLength
    # CDR2
    cdrLength = (cloneAnnot['cdr2.end'] - cloneAnnot['cdr2.start'] + 1) / 3
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['cdr2'] = cdrLength
    #CDR3
    cdrLength = (cloneAnnot['cdr3.end'] - cloneAnnot['cdr3.start'] + 1) / 3
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['cdr3'] = cdrLength
    #FR1
    cdrLength = (cloneAnnot['fr1.end'] - cloneAnnot['fr1.start'] + 1) / 3
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr1'] = cdrLength
    #FR2
    cdrLength = (cloneAnnot['fr2.end'] - cloneAnnot['fr2.start'] + 1) / 3
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr2'] = cdrLength
    #FR3
    cdrLength = (cloneAnnot['fr3.end'] - cloneAnnot['fr3.start'] + 1) / 3
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr3'] = cdrLength
    #FR4
    cdrLength = (cloneAnnot['fr4.end'] - cloneAnnot['fr4.start'] + 1) / 3
    cdrLength = Counter(cdrLength.tolist())
    spectraTypes['fr4'] = cdrLength
        
    return spectraTypes