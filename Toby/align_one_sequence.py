import sys
import os
import gzip
from seqfile import *
from dna import *
import string
import itertools
import functools
import collections
from plumb.bob import *
import cStringIO

def is_human((name, seq)):
    return name.split('|')[2] == 'Homo sapiens'

def is_human((name, seq)):
    return name.split('|')[2] == 'Homo sapiens'

def is_v_region((name, seq)):
    return name.split('|')[4] == 'V-REGION'

def is_functional((name, seq)):
    return name.split('|')[3] == 'F'

def is_type((name, seq), t):
    return name.split('|')[1].startswith(t)

DNA_SCORE = make_DNA_scoring_matrix(match=1, mismatch=-1, nmatch=0)

inf = cStringIO.StringIO(
'''\
>IGHV3-74*03_(L1)
GAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTAGTTCAGCCTGGGGGGTCCCTGAGACTC
TCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACTGGATGCACTGGGTCCGCCAAGCT
CCAGGGAAGGGGCTGGTGTGGGTCTCACGTATTAATAGTGATGGGAGTAGCACAACGTAC
GCGGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTAT
CTGCAAATGAACAGTCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCAAGAGA
''')

inf = open('Standard_IGHV_Repertoire.fa')
inf = fasta(inf)

inf = open('imgt.fa')
inf = list(itertools.ifilter(lambda x: is_type(x, 'IG') and is_human(x) and is_functional(x), fasta(inf)))

def assign(subject):
    out = []
    for name,query in inf:
        ali = local_align(
            subject, len(subject),
            query, len(query),
            DNA_MAP[0], DNA_MAP[1],
            DNA_SCORE,
            -5, -1)
        out.append((ali[0].score, name))
        alignment_free(ali)
    return max(out)

for name,seq in fasta(sys.stdin):
    r = assign(seq)

    print name, r[0], r[1].split('|')[1]
