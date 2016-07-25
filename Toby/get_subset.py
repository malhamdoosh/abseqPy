import sys
import os
import gzip
from seqfile import *
from dna import *
import string
import itertools
import functools
import collections
import cStringIO

def is_human((name, seq)):
    return name.split('|')[2] == 'Homo sapiens'

def is_human((name, seq)):
    return name.split('|')[2] == 'Homo sapiens'

def is_v_region((name, seq)):
    return name.split('|')[4] == 'J-REGION'

def is_functional((name, seq)):
    return name.split('|')[3] == 'F'

def is_type((name, seq), t):
    return name.split('|')[1].startswith(t)

inf = open('imgt.fa')
inf = list(itertools.ifilter(lambda x: is_type(x, 'IGL') and is_v_region(x) and is_human(x) and is_functional(x), fasta(inf)))

for a,b in inf:
  print '>{0}\n{1}'.format(a,b)
