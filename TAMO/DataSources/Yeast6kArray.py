"""
Yeast6kArray.py -- Information about the 6,000-probe intergenic array used in papers from the
                   laboratory of Richard Young.  Array uses PCR products, and much of the data
                   accessed by the routines herein are based on in silico PCR predictions by the
                   Fraenkel lab.

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
                   
"""

import re, math
from   TAMO.DataSources import SGD
import TAMO.paths

#SGD-like information, but not directly from SGD
_ORF2GENE  = {}
_GENE2ORF  = {}
_ORF2FUNC  = {}

#Simple 1-to-1 mapping.  (Misses, for example, divergent promoters)
_ORF2PROBE = {}
_PROBE2ORF = {}

#For more complex (but outdated) mapping.  Return lists
_ORF2PROBES = {}
_PROBE2ORFS = {}
_MAPTYPE    = {}


_SPOT_POS = {}
_SPOT_BIAS= {}
_PCR_BAD  = {}

_orf2probefile   = TAMO.paths.Yeast6kArraydir + 'orf2probe.csv'
_EFprobe2orffile = TAMO.paths.Yeast6kArraydir + 'probe2orfmapEF031109_uniqPCR.tsv'
_EForf2probefile = TAMO.paths.Yeast6kArraydir + 'orf2probemapEF031109_uniqPCR.tsv'
_probeposfile    = TAMO.paths.Yeast6kArraydir + 'YeastProbe05pcnt.csv'
_pcrfile         = TAMO.paths.Yeast6kArraydir + 'pcrresults.csv'
_biasfile        = TAMO.paths.Yeast6kArraydir + 'feature_bias_v7.csv'

TAMO.paths.CHECK([_orf2probefile, _EFprobe2orffile, _EForf2probefile, _probeposfile, _pcrfile, _biasfile],
                 'Whitehead')

def _loadfiles():
    global _orf2probefile
    global _ORF2PROBE, _PROBE2ORF, _ORF2GENE, _ORF2FUNC
    FID = open(_orf2probefile,'r'); lines = FID.readlines(); FID.close()
    for i in range(1,len(lines)):
        toks = lines[i].split(',')
        toks = [t.strip() for t in toks]
        orf, probe, gene, function = toks[0], toks[2], toks[8], re.sub('"','',toks[9])
        _ORF2PROBE[orf]   = probe
        _PROBE2ORF[probe] = orf
        _ORF2GENE[orf]    = gene
        _GENE2ORF[gene]   = orf
        _ORF2FUNC[orf]    = function
    #_load_Xfiles()
    _load_EFfiles()

def _load_EFfiles():
    FID = open(_EForf2probefile,'r'); o2plines = FID.readlines(); FID.close()
    FID = open(_EFprobe2orffile,'r'); p2olines = FID.readlines(); FID.close()
    for line in p2olines:
        toks  = [x.strip() for x in line.split()]
        probe = toks[0]
        itype = toks[1]
        orfs = [re.sub('[AB]\|','',x) for x in toks[2:]]
        _PROBE2ORFS[probe] = orfs
    for line in o2plines:
        toks  = [x.strip() for x in line.split('\t')]
        orf   = toks[0]
        otype = toks[1]
        odist = int(toks[2])
        itype = toks[3]
        probes= [re.sub('[AB]\|','',x.split()[0]) for x in toks[4:] if re.search('[AB]\|',x)]
        _ORF2PROBES[orf] = probes
        
                
                   

def _load_Xfiles():
    FID = open(_Xorf2probefile,'r'); lines = FID.readlines(); FID.close()
    for line in lines[1:]:
        toks = line.split(',')
        orf, probe, category = toks[0], toks[1], re.sub('"','',toks[7])
        if not _ORF2PROBES.has_key(orf): _ORF2PROBES[orf] = []
        _ORF2PROBES[orf].append(probe)
        if not _PROBE2ORFS.has_key(probe): _PROBE2ORFS[probe] = []
        _PROBE2ORFS[probe].append(orf)
        if not _MAPTYPE.has_key(orf): _MAPTYPE[orf] = {}
        _MAPTYPE[orf][probe] = category

def _result(D,key,fail=None):
    if not D: _loadfiles()
    if not D.has_key(key):
        if fail==None: fail=key
        return fail
    else:                    return D[key]

def orf2gene(orf):     return _result(_ORF2GENE, orf)
def orf2probe(orf):    return _result(_ORF2PROBE,orf)
def probe2orf(probe):  return _result(_PROBE2ORF,probe)
def gene2orf(gene):    return _result(_GENE2ORF, gene)
def orf2func(orf):     return _result(_ORF2FUNC, orf)

def probe2orfs(probe): return _result(_PROBE2ORFS,probe,[])
def orf2probes(orf):   return _result(_ORF2PROBES,orf,[])
def probe2genes(probe):return [orf2gene(x) for x in probe2orfs(probe)]
def probe2funcs(probe):return [orf2func(x) for x in probe2orfs(probe)]
def gene2probes(gene): return orf2probes(gene2orf(gene))

def probe2func(probe): return orf2func ( probe2orf(probe))
def gene2func(gene)  : return orf2func (  gene2orf(gene ))
def gene2probe(gene) : return orf2probe(  gene2orf(gene ))
def probe2gene(probe): return orf2gene ( probe2orf(probe))

def maptype(orf,probe):
    if not _MAPTYPE: _loadfiles()
    if _MAPTYPE.has_key(orf):
        if _MAPTYPE[orf].has_key(probe):
            return _MAPTYPE[orf][probe]
        else:
            return ''
    elif _MAPTYPE.has_key(probe):
        if _MAPTYPE[probe].has_key(orf):
            return _MAPTYPE[probe][orf]
        else:
            return ''
    return ''


############################### STUFF for SPOTS on the ARRAY #######################
class chr_pos:
    """
    A chromosomal position
    """
    def __init__(self,start=None,stop=None,chr=None):
        self.start = start
        self.stop  = stop
        self.chr   = chr
    def __repr__(self):
        s = '%d : %8d - %8d'%(self.chr,self.start,self.stop)
        return s

class Spot:
    """
    A spot on the microarrray (name, position(s))
    """    
    def __init__(self,csvline=''):
        self.name      = None
        self.npositions= 0
        self.positions = []
        if csvline:
            toks = csvline.split(',')
            self.name = toks[1]
            self.npositions = int(toks[8])
            if int(toks[8]) > 0:
                locs = toks[9].split('][')
                for loc in locs:
                    toks = re.sub('[\[\]a-z\:\(\)\-]',' ',loc).split()
                    chr, start, stop = int(toks[0]), int(toks[1]), int(toks[2])
                    self.positions.append(chr_pos(start,stop,chr))
                    
def spot_pos(probe_id):
    """
    What are the chromosomal coordinates associated with the ePCR prediction of this spot?
    """
    if not _SPOT_POS:
        F = open(_probeposfile,'r'); lines = F.readlines(); F.close()
        del lines[0] #Header lines
        for line in lines:
            _spot = Spot(line.strip())
            _SPOT_POS[_spot.name] = _spot
    if _SPOT_POS.has_key(probe_id):
        return _SPOT_POS[probe_id].positions
    else:
        return []
    
def spot_bias(probe_id):
    if not _SPOT_BIAS:
        F = open(_biasfile,'r'); lines = F.readlines(); F.close()
        for line in lines:
            toks = line.strip().split(',')
            if len(toks) < 2: continue
            probe, biastxt = toks
            if not biastxt: continue
            _SPOT_BIAS[probe] = float(biastxt)
    if _SPOT_BIAS.has_key(probe_id):
        return _SPOT_BIAS[probe_id]
    else:
        return 0

def pcr_bad(probe_id):
    if not _PCR_BAD:
        F = open(_pcrfile,'r'); lines = F.readlines(); F.close()
        del lines[0]
        for line in lines:
            toks = line.split(',')
            probe = toks[3]
            badness  = ''
            goodness = ''
            for testdate in toks[5:]:
                subtoks = testdate.split()
                if len(subtoks) > 1:
                    badness = badness + subtoks[1]
                elif len(subtoks) < 1:
                    badness = badness + '-'
                elif subtoks[0] == 'n/a':
                    badness = badness + 'N'
                else:
                    goodness = 'Y'
            if goodness:
                _PCR_BAD[probe] = ''
            else:
                _PCR_BAD[probe] = badness
    if not _PCR_BAD.has_key(probe_id):
        return None
    else:
        return _PCR_BAD[probe_id]
        
def telomeredist(spotname):
    """
    Distance from feature (spot) to nearest telomere
    """
    poss = spot_pos(spotname)
    if poss:
        teldist =  min([SGD.dist_from_tel(pos.chr,pos.start,pos.stop) for pos in poss])
        return teldist
    return -1

def cendist(spotname):
    """
    Distance from feature (spot) to nearest centromeme
    """
    poss = spot_pos(spotname)
    if poss:
        cendist =  min([SGD.dist_from_cen(pos.chr,pos.start,pos.stop) for pos in poss])
        return cendist
    return -1

def stringent_filter(ids,bias_cutoff=10):
    result = []
    for id in ids:
        if (pcr_bad(id)             or
            #math.fabs(spot_bias(id))  >  bias_cutoff  or
            #telomeredist(id)  <  100 or
            #cendist(id)       <  100 or
            len(spot_pos(id)) > 2):
            continue
        result.append(id)
    return result
