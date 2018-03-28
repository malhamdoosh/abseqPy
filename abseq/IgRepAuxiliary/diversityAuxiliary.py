'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from collections import Counter, defaultdict
from numpy import logical_not

_REGIONS = ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4']


#  spectratype, that is, histogram of clone counts by CDR/FR nucleotide length.
#  The spectratype is useful to detect pathological and highly clonal repertoires, 
#  as the spectratype of non-expanded T- and B-cells has a symmetric gaussian-like distribution.
def annotateSpectratypes(cloneAnnot, amino=True):
    # TODO: add annotation to clonotypes, e.g., germline genes
    # TODO: add nucleotide level calculations
    denom = 3 if amino else 1
    spectraTypes = {}
    for region in _REGIONS:
        spectraType = ((cloneAnnot[region + '.end'] - cloneAnnot[region + '.start'] + 1) / denom).astype(int)
        spectraTypes[region] = Counter(spectraType)
    # V domain
    spectraType = ((cloneAnnot['fr4.end'] - cloneAnnot['fr1.start'] + 1) / denom).astype(int)
    spectraTypes['v'] = Counter(spectraType)

    return spectraTypes


# clonotype is the histogram of clone counts by CDR/FR amino acid sequence, partitioned by V germline (gene level)
def annotateClonotypes(cloneSeqs, removeNone=True):
    clonoTypes = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for row in cloneSeqs.itertuples():
        geneName = row.germline.split("*")[0]
        variableAA = ""
        badAA = False
        for region in _REGIONS:
            regionAA = getattr(row, region)
            variableAA += regionAA
            if regionAA == 'None':
                badAA = True
            if not (removeNone and regionAA == "None"):
                clonoTypes[geneName][region][regionAA] += 1
        if not (badAA and removeNone):
            clonoTypes[geneName]['v'][variableAA] += 1

    return clonoTypes
