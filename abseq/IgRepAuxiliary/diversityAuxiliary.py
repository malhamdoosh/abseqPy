'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
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
    spectraType = ((cloneAnnot['cdr1.end'] - cloneAnnot['cdr1.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['cdr1'] = spectraType
    # CDR2
    spectraType = ((cloneAnnot['cdr2.end'] - cloneAnnot['cdr2.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['cdr2'] = spectraType
    # CDR3
    spectraType = ((cloneAnnot['cdr3.end'] - cloneAnnot['cdr3.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['cdr3'] = spectraType
    # FR1
    spectraType = ((cloneAnnot['fr1.end'] - cloneAnnot['fr1.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['fr1'] = spectraType
    # FR2
    spectraType = ((cloneAnnot['fr2.end'] - cloneAnnot['fr2.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['fr2'] = spectraType
    # FR3
    spectraType = ((cloneAnnot['fr3.end'] - cloneAnnot['fr3.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['fr3'] = spectraType
    # FR4
    spectraType = ((cloneAnnot['fr4.end'] - cloneAnnot['fr4.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['fr4'] = spectraType
    # V domain
    spectraType = ((cloneAnnot['fr4.end'] - cloneAnnot['fr1.start'] + 1) / 3).astype(int)
    spectraType = Counter(spectraType.tolist())
    spectraTypes['v'] = spectraType
        
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
    # V domain
    seqs = map(lambda x:''.join(x), zip(cloneSeqs['fr1'].tolist(),
                                        cloneSeqs['cdr1'].tolist(),
                                        cloneSeqs['fr2'].tolist(),
                                        cloneSeqs['cdr2'].tolist(),
                                        cloneSeqs['fr3'].tolist(),
                                        cloneSeqs['cdr3'].tolist(),
                                        cloneSeqs['fr4'].tolist()                                        
                                        )
               )
    clonoTypes['v'] = Counter(seqs)
    return clonoTypes











