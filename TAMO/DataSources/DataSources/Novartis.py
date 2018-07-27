"""
Novartis.py

Loads Novartis mapping Data and Dataset as a HT.Dataset object

id2feature(id)
id2other(id,other), where "other" can be 

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon

"""
from   gzip        import GzipFile
from   TAMO.HT     import Dataset
import TAMO.paths

ANNOFILE = TAMO.paths.Novartisdir + 'gnf1h-anntable.txt.gz'
DATAFILE = TAMO.paths.Novartisdir + 'U133A+GNF1B_101402.AD.txt'
ANNO     = {}
REVANNO  = {}
# feature    name    location     EC     desc
# 0          1        2       3          4                                  5                                6    7                         8
#'1005_at', 'DUSP1', '5q34', '3.1.3.16', 'dual specificity phosphatase 1', 'dual specificity phosphatase 1', '', 'Cell stress;Hydrolase;', 'oxidative stress response;non-membrane spanning protein tyrosine
#                RefSeq                  LocusLink SwissProt Unigene
#                 9           10  11  12  13      14        15
# phosphatase;', 'NM_004417', '', '', '', '1843', 'P28562', 'Hs.171695'

#0              1       2               3                                               4                                                               5                               6               7       8
#Taxon	        Name	Probeset ID	Reporters	                                Genome Location	                                                LocusLink RefSeq	        UniGene	        UniProt	Ensembl
#Homo sapiens	XG	gnf1h08239_s_at	gnf1h02393_at;gnf1h08238_at;gnf1h08239_s_at	ChrX:2.265-2.330 (+) (NCBI34);ChrY:2.265-2.284 (+) (NCBI34)	7499	NM_175569;NP_780778	Hs.179675	P55808	ENSP00000262688

#9      10                                                                      11                                                                                                              12   
#Aliase	Description	                                                        Function	                                                                                                Protein Families
#PBDX	Xg blood group (pseudoautosomal boundary-divided on the X chromosome)	DNA binding (GO:0003677);biological_process unknown (GO:0000004);molecular_function unknown (GO:0005554)	Proline-rich region (IPR000694)


ref = {'feature':2,
       'name':1,
       'location':4,
       'desc':10,
       'func':11,
       'refseq':5,
       'locuslink':5,
       'uniprot':7,
       'unigene':6,
       'ensenbl':8,
       'family':12}

def load_anno():
    if ANNO: return
    TAMO.paths.CHECK(ANNOFILE,'Novartis')
    F     = GzipFile(ANNOFILE)
    lines = F.readlines()
    F.close()
    for line in lines:
        toks = line.strip().split('\t')
        if len(toks) < 15: continue
        D={}
        for id,pos in ref.items():
            id = id.lower()
            if pos != 0: REVANNO[toks[pos]] = toks[0]
            D[id] = toks[pos]
            ANNO[toks[0]] = D

def id2feature(id):
    load_anno()
    if REVANNO.has_key(id):
        return REVANNO[id]
    else:
        return ''

def id2other(id,other):
    other = other.lower()
    feature = id2feature(id)
    if not feature: return ''
    if not ANNO[feature].has_key(other): return ''
    return ANNO[feature][other]

def humandata():
    TAMO.paths.CHECK(DATAFILE,'Novartis')
    G = Dataset(DATAFILE)

def id2unigene(id):   return id2other(id,'unigene')
def id2locuslink(id): return id2other(id,'locuslink')
def id2ll(id):        return id2locuslink(id)
def id2uniprot(id):   return id2other(id,'uniprot')

                



