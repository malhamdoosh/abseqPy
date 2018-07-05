'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from collections import Counter, defaultdict

_REGIONS = ['fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4']


#  spectratype, that is, histogram of clone counts by CDR/FR nucleotide length.
#  The spectratype is useful to detect pathological and highly clonal repertoires, 
#  as the spectratype of non-expanded T- and B-cells has a symmetric gaussian-like distribution.
def annotateSpectratypes(cloneAnnot, amino=True):
    # TODO: add annotation to clonotypes, e.g., germline genes
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
def annotateClonotypes(cloneSeqs, segregate=False, removeNone=True):
    if segregate:
        clonoTypes = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

        for row in cloneSeqs.itertuples():
            geneName = row.germline.split("*")[0]
            variableAA = ""
            badAA = False
            for region in _REGIONS:
                regionAA = getattr(row, region)
                variableAA += regionAA
                if regionAA == 'None' or regionAA == '':
                    badAA = True
                if not (removeNone and (regionAA == "None" or regionAA == "")):
                    clonoTypes[geneName][region][regionAA] += 1
            if not (badAA and removeNone):
                clonoTypes[geneName]['v'][variableAA] += 1
    else:
        clonoTypes = {}

        for region in _REGIONS:
            seqs = cloneSeqs[region].tolist()
            clonoTypes[region] = Counter(seqs)
            if removeNone:
                clonoTypes[region].pop("None", None)
                clonoTypes[region].pop("", None)

        def _join(frcdr):
            if removeNone and ('None' in frcdr or '' in frcdr):
                return "None"
            return "".join(frcdr)

        # V domain
        seqs = map(_join, zip(cloneSeqs['fr1'].tolist(),
                              cloneSeqs['cdr1'].tolist(),
                              cloneSeqs['fr2'].tolist(),
                              cloneSeqs['cdr2'].tolist(),
                              cloneSeqs['fr3'].tolist(),
                              cloneSeqs['cdr3'].tolist(),
                              cloneSeqs['fr4'].tolist()
                              )
                   )
        clonoTypes['v'] = Counter(seqs)
        if removeNone:
            clonoTypes['v'].pop("None", None)
            clonoTypes['v'].pop("", None)

    return clonoTypes

