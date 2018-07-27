#!env python -O
"""
GO.py -- Access and statistis for yeast GO-slim

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
"""
import sys, re, os, math, time, pickle

import TAMO.paths
from   TAMO.DataSources import SGD
from   TAMO.util        import Arith

_GO_ALL_FILE = TAMO.paths.SGDdir + 'go_slim_mapping.tab'
TAMO.paths.CHECK(_GO_ALL_FILE,'GO')
#ftp://genome-ftp.stanford.edu/pub/yeast/data_download/literature_curation/

BAD_CATEGORIES = []

class Annotation:
    """
    Data structure to store a GO-slim entry
    """
    def __init__(self,go_line=None):
        if go_line: self._parse_line(go_line)
    def _parse_line(self,line):
        toks = [x.strip() for x in line.split('\t')]
        self.orf      = toks[0]
        self.aspect   = toks[3]
        self.desc     = toks[4]
        self.isnot    = 0   #(toks[6] == 'Y')  #Commented to accomodate go_slim_mapping
        self.evidence = ''  #toks[7]
        self.pubmed   = ''  #toks[9]
    def __repr__(self):
        s = '<GO: %s %s (%s %s)>'%(
            self.aspect, self.desc, self.evidence, self.pubmed)
        return s

        
#Global Variables

_gos_by_orf   = {}   #Dictionary (by orf) of lists of annotations
_orfs_by_go   = {}   #Dictionary (by go description) of lists of orf names
_by_size      = {}   #Dictionary (by number of genes) describing distribution of category sizes

def load_GO():
    print "#Loading GO"
    FID = open(_GO_ALL_FILE,'r'); lines = FID.readlines(); FID.close()
    for line in lines:
        anno = Annotation(line)
        #if anno.aspect == 'C': continue  #Ignore all but "Process"
        if anno.aspect != 'P':             continue  #Ignore all but "Process"
        if anno.desc    in BAD_CATEGORIES: continue
        if anno.isnot:                     continue #Ignore "negative" GO ontologies
        if not _orfs_by_go.has_key(anno.desc): _orfs_by_go[anno.desc] = []
        _orfs_by_go[anno.desc].append(anno.orf)
        if not _gos_by_orf.has_key(anno.orf): _gos_by_orf[anno.orf] = []
        _gos_by_orf[anno.orf].append(anno)
    for orflist in _orfs_by_go.values():
        size = len(orflist)
        if not _by_size.has_key(size): _by_size[size] = 0.0
        _by_size[size] = _by_size[size] + 1.0

# List annotations for a gene
def gene2annotations(gene):  return annotations(SGD.gene2orf(gene))

# List annotations for an ORF
def annotations(orf):
    if not _gos_by_orf: load_GO()
    if not _gos_by_orf.has_key(orf):
        return []
    return _gos_by_orf[orf]

# List all orfs associated with an annotation
def annotation2orfs(anno_or_txt):
    if not _orfs_by_go: load_GO()
    if type(anno_or_txt) == type(''):
        key = anno_or_txt
    else:
        key = anno_or_txt.desc
    if not _orfs_by_go.has_key(key):
        return []
    else:
        return _orfs_by_go[key]

# Which annotations for a gene are also assocaited with probes bound in the
# Yeast6k array (e.g. from a particular experiment)
def geneprobelistmatch(gene,probelist):
    if 'Yeast6kArray' not in dir(): from TAMO.DataSources import Yeast6kArray
    s_orf = SGD.gene2orf(gene)
    orfs = []
    for probe in probelist:
        orfs.extend(Yeast6kArray.probe2orfs(probe))
    #print "%4d --> %4d   "%(len(probelist),len(orfs)),
    return orforflistmatch(s_orf,orfs)
    
def orforflistmatch(s_orf,orflist):
    """
    Which annotations associated with s_orf are significantly overrepresented in orflist?
    """
    if not _orfs_by_go: load_GO()
    norfs   = len(orflist)
    totorfs = len(_gos_by_orf.keys())

    categories = []
    for anno in annotations(s_orf):
        categories.append(anno.desc)

    sigs = []
    for category in categories:
        all_by_cat = annotation2orfs(category)
        nall       = len(all_by_cat)
        fracall    = float(nall)/float(totorfs)
        
        sub_by_cat = [x for x in orflist if x in all_by_cat]
        nsub       = len(sub_by_cat)

        if norfs == 0: fracsub = 0
        else:          fracsub    = float(nsub)/float(norfs)

        sig        = Arith.hypgeomsummore(nall,totorfs,norfs,nsub)
        sigs.append( (sig, category, fracall, fracsub, sub_by_cat) )
        #if s_orf == SGD.gene2orf('NRG1'):
        #    print '@ NRG1 %5.2e %-40s %4d %4d %4d %4d '%(
        #        sig,category,nall,totorfs,norfs,nsub)
    sigs.sort()

    #for sig in sigs: print sig
    return sigs

def probelist2categories(probelist,thresh=0.05):
    """
    Which categories are overrepresented among probes bound in in a Yeast6k array experiment?
    """
    if 'Yeast6kArray' not in dir(): from TAMO.DataSources import Yeast6kArray
    orfs = []
    preorfs = []
    for probe in probelist:
        preorfs.extend(Yeast6kArray.probe2orfs(probe))
    for o in preorfs:
        if (not o in orfs): orfs.append(o)
    return orflist2categories(orfs,thresh)

def orfs2cats(orflist,thresh=0.05):
    """
    [[Same as orflist2categories]]

    Which categories are overrepresented among orflist.  Thresh is applied after Bonferroni
    correction.

    Returns: sorted list of tuples, [(signifance, category), (significance, category), ...]
    """
    return orflist2categories(orflist,thresh)

def orflist2categories(orflist,thresh=0.05):
    """
    Which categories are overrepresented among orflist.  Thresh is applied after Bonferroni
    correction.

    Returns: sorted list of tuples, [(signifance, category), (significance, category), ...]
    """
    ans = []
    for tup in orflist2categories_long(orflist, thresh):
        ans.append( (tup[0],tup[1]))
    return ans

def orflist2categories_long(orflist,thresh=0.05):
    """
    Which categories are overrepresented among orflist.  Thresh is applied after Bonferroni
    correction.

    Returns: Sorted list of tuples, [(signifance, category, fracall, fracsub, orflist), ...]
    """
    if not _orfs_by_go: load_GO()
    norfs   = len(orflist)
    totorfs = len(_gos_by_orf.keys())

    #Determine which categories are described by the set
    categories = []
    for orf in orflist:
        for anno in annotations(orf):
            if anno.desc not in categories:
                categories.append(anno.desc)

    totcats = float(len(categories))
    sigs = []
    for category in categories:
        all_by_cat = annotation2orfs(category)
        nall       = len(all_by_cat)
        fracall    = float(nall)/float(totorfs)
        
        sub_by_cat = [x for x in orflist if x in all_by_cat]
        nsub       = len(sub_by_cat)
        fracsub    = float(nsub)/float(norfs)

        sig        = Arith.hypgeomsummore(nall,totorfs,norfs,nsub) * totcats
        #print category,sig,nall,totorfs,norfs,nsub,'       ',totcats
        sigs.append( (sig, category, fracall, fracsub, sub_by_cat) )
        
    sigs.sort()

    ans = []
    for sigdata in sigs:
        sig, category, fracall, fracsub, orfs = sigdata
        if sig > thresh: continue
        ans.append(sigdata)

    return ans


