"""
SGD.py -- Routines for convenient access to SGD data, including cevivsiae genome sequence.


Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
"""
import sys, re, math, pickle
import TAMO.paths
from TAMO.seq   import Fasta
from TAMO.paths import SGDdir

_ORFFILE  = SGDdir + 'SGD_features.tab'
_EXTID    = SGDdir + 'dbxref.tab'
_SWPFASTA = SGDdir + 'yeast_nrpep.fasta.gz'
_CHR_LEN  = SGDdir + 'chromosome_length.tab'
_LOCALIZE = SGDdir + 'Huh_Nature_2003.tab'
_ORFPSEQS = SGDdir + 'orf_trans_all.fasta.gz'

TAMO.paths.CHECK([_ORFFILE, _EXTID, _SWPFASTA, _CHR_LEN, _LOCALIZE, _ORFPSEQS], 'SGD')

_featuresD = {}
_gene2orf  = {}
_orf2gene  = {}
_orf2swp   = {}
_orf2ext   = {}
_swp_seqs  = {}
_swp2swp   = {}
_localize  = {}
_abundance = {}
_orfpseqs  = {}



'''
The dbxref file has entries like this:

CA3405      CandidaDB Gene ID                    YJL221C S000003757
1787        DIP       Gene ID                    YJL221C S000003757
P40439      SIB       Swiss-Prot ID              YJL221C S000003757
853235      NCBI      Gene ID                    YJL221C S000003757
Z34098      NCBI      DNA version ID             YJL221C S000003757
AB109221    NCBI      DNA version ID             YJL221C S000003757
Z49496      NCBI      DNA version ID             YJL221C S000003757
D43761      NCBI      DNA version ID             YJL221C S000003757
BAD00094.1  NCBI      Protein version ID         YJL221C S000003757
CAA83990.1  NCBI      Protein version ID         YJL221C S000003757
CAA89517.1  NCBI      Protein version ID         YJL221C S000003757
BAA07818.1  NCBI      Protein version ID         YJL221C S000003757
NP_012314.1 NCBI      RefSeq protein version ID  YJL221C S000003757
3.2.1.20    IUBMB     EC number                  YJL221C S000003757

'''

def orf2other(qorf,qtype):
    global _EXTID
    global _orf2ext
    if not _orf2ext:
        FID = open(_EXTID,'r')
        lines = FID.readlines()
        FID.close()
        for line in lines:
            toks = [x.strip() for x in line.split('\t')]
            if toks[1] in ['SIB', 'DIP', 'IUBMB', 'CandidaDB']:
                exid,type,orf = toks[0], toks[1], toks[-2]
                if type == 'SIB': type = 'SwissProt'
            elif (toks[1] == 'NCBI') and (toks[2] == 'RefSeq protein version ID'):
                exid,type,orf = toks[0],'RefSeq',toks[-2]
            else: continue
            try:
                _orf2ext[orf][type] = exid
            except:
                _orf2ext[orf] = {}
                _orf2ext[orf][type] = exid
                    
    try:
        return(_orf2ext[qorf][qtype])
    except:
        return(None)
                
def get_swp_sequence(query):
    global special_cases
    orf, swp = None, None
    if special_cases.has_key(query):
        orf = gene2orf(special_cases[query])
    else:
        orf = gene2orf(query)
    if orf:
        swp = orf2swp(orf)
    else:
        swp = orf2swp(query)
    if not swp:
        print '# Could not find SwissProt ID for %s (%s)'%(query,orf)
        return None
    txt = swp_find_and_format(swp)
    if not txt:
        print '# Could not find Sequence for Swissprot ID %s (query %s)'%(swp,query)
        return None
    ans = '>%-6s %-7s %-12s %-12s\n%s'%(query,orf,swp2swp(swp),swp,txt)
    return (ans,orf)

def swp_find_and_format(swp):
    global _swp_seqs
    if not _swp_seqs:
        _swp_seqs = Fasta.load(_SWPFASTA,key_func=lambda x:x)
    hits = []
    for key in _swp_seqs.keys():
        if key[0:60].find(swp) >= 0:
            hits.append(key)
    if not hits:
        return None
    if len(hits) > 1:
        print "# Multiple matches found for %s:"%swp
        for hit in hits: print '#',hit
        return None
    hit = hits[0]
    seq = _swp_seqs[hit]
    txt = ''
    for i in range(0,len(seq),70):
        txt = txt + seq[i:i+70] + '\n'
    return txt
    

def experiment2gene(experiment):
    global special_cases
    gene = re.split('[ _]*',experiment)[0]
    if special_cases.has_key(gene):
        gene = special_cases[gene]
    return gene

def experiment2orf(experiment):
    return gene2orf(experiment2gene(experiment))

def swp2swp(swp):
    'Converts, when possible, from P014543 to ADR1_YEAST'
    'Only works for yeast right now'
    global _swp2swp
    if not _swp2swp:
        lines = Fasta.keys(_SWPFASTA,key_func=lambda x:x)
        for line in lines:
            toks = line.split()
            text_name  = toks[1]
            numeric_name = toks[2]
            if text_name[0:2] == 'SW' and numeric_name[0] == 'P':
                _swp2swp[text_name[3:]]  = numeric_name
                _swp2swp[numeric_name]   = text_name[3:]
    if _swp2swp.has_key(swp):
        return _swp2swp[swp]


def gene2swp(gene):   return(swp2swp(orf2swp(gene2orf(gene))))

def experiment2swp(experiment):
    gene = experiment2gene(experiment)
    return swp2swp(orf2swp(gene2orf(gene)))

def orf2swp(orf_query):
    if not _orf2swp:
        FID = open(_ORFSWP,'r')
        lines = FID.readlines()
        FID.close()
        for line in lines:
            toks = line.split('\t')
            swp,orf = toks[0],toks[2]
            _orf2swp[orf] = swp
    if _orf2swp.has_key(orf_query):
        return _orf2swp[orf_query]
    else:
        return None

def orf2swp(orf_query):
    return orf2other(orf_query,'SwissProt')



_chr_lengths = {}
def chr_len(number):
    if not _chr_lengths:
        FID = open(_CHR_LEN,'r'); lines = FID.readlines(); FID.close()
        for line in lines:
            toks = [x.strip() for x in line.split('\t')]
            chr, length = int(toks[0]), int(toks[2])
            _chr_lengths[chr] = length
    if _chr_lengths.has_key(int(number)):
        return _chr_lengths[int(number)]
    else:
        return None


def feature_dist(Astart,Astop, Bstart, Bstop):
    FABS = math.fabs
    dist = min( FABS(Astart - Bstart), FABS(Astart-Bstop),
                FABS(Astop  - Bstart), FABS(Astop-Bstop) )
    return dist

def feature_pos(_feature):
    F = feature(_feature)
    if not F: F = feature(gene2orf(_feature))
    if not F: return None
    else:     return (F.chr,F.start,F.stop)

def dist_from_cen(chr,start,stop):
    if chr < 10: censtr = 'CEN%d'%chr
    else       : censtr = 'CEN%d'%chr
    cen = feature(censtr)
    return feature_dist(cen.start, cen.stop, start, stop)

def dist_from_tel(chr,start,stop):
    return min(start, chr_len(chr) - stop)

def _load_features():
    global _gene2orf, _orf2gene
    global _featuresD
    FID = open(_ORFFILE,'r')
    lines = FID.readlines()
    FID.close()
    for line in lines:
        toks = [x.strip() for x in line.split('\t')]
        if toks[1][0:3] == 'ORF' and toks[1]:
            gene,orf = toks[4],toks[3]
            _gene2orf[gene] = orf
            _orf2gene[orf]  = gene
        try:
            name = toks[3]
            alts = toks[5].split('|')
            ftype= toks[1]
            chr  = int(toks[8])
            start= int(toks[9])
            stop = int(toks[10])
            desc = toks[15]
            entry = featureline(name,alts,ftype,chr,start,stop,desc)
            if not _featuresD.has_key(ftype): _featuresD[ftype] = {}
            _featuresD[ftype][name] = entry
        except: continue
            
            
class featureline:
    def __init__(self,name,alts,ftype,chr,start,stop,desc):
        self.name = name
        self.alts = alts
        self.type = ftype
        self.chr  = chr
        self.start= start
        self.stop = stop
        self.desc = desc
    def __repr__(self):
        s = '%-10s %s %s %2d %8d %8d : %s'%(
            self.name, '|'.join(self.alts), self.type, self.chr, self.start, self.stop, self.desc)
        return s

def gene2orf(gene_query):
    global special_cases
    global _gene2orf
    if not _gene2orf:
        _load_features() 
    if _gene2orf.has_key(gene_query):
        return _gene2orf[gene_query]
    elif (gene_query in special_cases.keys()):
        if _gene2orf.has_key(special_cases[gene_query]):
            return _gene2orf[special_cases[gene_query]]
        else:
            return special_cases[gene_query]
    elif (gene_query in _gene2orf.values()):
        return gene_query
    else:
        return gene_query
    
def orf2gene(orf_query):
    if not _orf2gene:
        _load_features()
    if _orf2gene.has_key(orf_query):
        return _orf2gene[orf_query]
    elif (orf_query in _orf2gene.values()):
        return orf_query
    else:
        return None
    
def feature(featurename):
    if not _featuresD: _load_features()
    for key,D in _featuresD.items():
        if D.has_key(featurename):
            return D[featurename]
    return None

def featuredict(featuretype):
    if not _featuresD: _load_features()
    if _featuresD.has_key(featuretype):
        return _featuresD[featuretype]
    return none

def gene2pseq(gene): return orf2pseq(gene2orf(gene))

def orf2pseq(orf):
    global _orfpseqs
    if not _orfpseqs:
        from TAMO.seq import Fasta
        _orfpseqs = Fasta.load(_ORFPSEQS)
        for _orf,pseq in _orfpseqs.items():
            if pseq[-1] == '*': _orfpseqs[_orf] = pseq[:-1]
    if _orfpseqs.has_key(orf): return _orfpseqs[orf]
    else:                      return ''



#LOCALIZATION/Abundance Routines

def genenonnuclear(gene):   return orfnonnuclear(gene2orf(gene))
def genelocalization(gene): return orflocalization(gene2orf(gene))
def geneabundance(gene): return orfabundance(gene2orf(gene))

def orfabundance(orf):
    if not _abundance: load_localization()
    if not _abundance.has_key(orf): return None
    return _abundance[orf]


def orfnonnuclear(orf):
    if not _localize: load_localization()
    if not _localize.has_key(orf):  return(None)
    nuc = 0 
    for loc in _localize[orf]:
        if loc.find('nuc') >= 0 and loc.find('non')<0 :
            nuc = 1
    return not nuc
    

def orflocalization(orf):
    if not _localize: load_localization()
    if not _localize.has_key(orf):  return(None)
    U = {}
    for loc in _localize[orf]:
        if not U.has_key(loc): U[loc]=0
        U[loc] = U[loc]+1
    #txt = '/'.join(['%s_%d'%(l,c) for l,c in U.items()])
    txt = ','.join(U.keys())
    return txt

def load_localization_Kumar():  #Not used as of release of Huh...O'shea data 10-17-03
    if not _localize:
        FID = open(_LOCALIZE,'r'); lines = FID.readlines(); FID.close()
        for line in lines[1:]:
            D={}
            toks = line.split(',')
            orf  = toks[0].strip()
            D['locs'] = [x.strip() for x in toks[2].split('*')]
            D['com']  = [x.strip() for x in toks[3].split('*')]
            _localize[orf] = D

def load_localization():
    if not _localize:
        FID = open(_LOCALIZE,'r'); lines = FID.readlines(); FID.close()
        labels = [x.strip() for x in lines[0].split('\t')]
        for line in lines[1:]:
            toks = [x.strip() for x in line.split('\t')]
            orf, name, abund, locs = toks[1], toks[2], toks[6], toks[8]
            _localize[orf]  = locs.split(',')
            if abund.isdigit(): _abundance[orf] = int(abund)
            else:               _abundance[orf] = -1


ChrD = {}
def get_seq(chr,start=None,stop=None):
    global ChrD
    if not ChrD:
        from TAMO.seq import Fasta
        ChrD = Fasta.load(SGDdir + 'NCBI_yeast_genome.fsa')
    if (type(chr) != type('')) or (chr.find('chr') != 0):  # 1 -> chr1, 'X' -> chrX
        chr = 'chr%s'%chr
    if (start == None) and chr.find(':') > 0:                  # chr4:454-465 -> chr4, 454, 465
        _chr,_range = chr.split(':')
        chr = _chr
        start, end = _range.split('-')
        start, end = int(start), int(end)
    return ChrD[chr][start-1:end]

special_cases = {'GRF10(Pho2)':'PHO2',
                 'SIG1': 'MOT2',
                 'YAP8': 'ARR1',
                 'HDT80':'NDT80',
                 'MATa1':'HMRA1',
                 'MATA1': 'HMRA1',
                 'BYE1': 'YKL005C',
                 'MIG3': 'YER028C',
                 'NNF2': 'YGR089W',
                 'USV1':'YPL230W'}

#    'MSN2d4':None
#    'MSN4d2':None
    



