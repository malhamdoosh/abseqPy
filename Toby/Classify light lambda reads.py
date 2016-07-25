import sys
import os
import gzip
from seqfile import *
from dna import *
import string
import itertools
from plumb.bob import *

def fq_hist(seqs):
    x = collections.Counter(seqs)
    y = collections.Counter(x.values())

    print 'sequence frequency histogram'
    for a,b in sorted(y.items()):
        print '{0:5d} {1:10d}'.format(a, b)
    

framework = [
  'QYELTQPPSASGTPGQRVTFSC',
  'WYQQLPGTAPKLLIY',
  'GVPDRFSASKSGTSASLAISGLQSEDEADYYC',
  'FGGGTQ'
]

def best_score(seq, query, min_pos = 0):
  seq_ = seq[min_pos:]
  alignment = glocal_align(seq_, len(seq_),
                           query, len(query),
                           PROT_MAP[0], PROT_MAP[1],
                           BLOSUM62,
                           -100, -100)

  hsps = []
  
  frag = alignment[0].align_frag
  while frag:
    frag = frag[0]
    if frag.type == MATCH:
      hsps.append(frag)
    frag = frag.next

  score = alignment[0].score
  match_pos = min_pos + hsps[0].sa_start
  match_len = hsps[-1].sa_start + hsps[-1].hsp_len - hsps[0].sa_start

  alignment_free(alignment)

  return score, match_len, match_pos

def is_light(s):
  return s[1][:6] in ('AAGACA', 'AGAGCG', 'AGTACG', 'AGGACA')

def first(n, iter):
  for i in xrange(n): yield iter.next()

N = 0
for name,seq,qual in itertools.ifilter(is_light, fastq(gzip.open('03_SeqPrep_merged/phage_merged.fq'))):
    seq = 'C'+seq

    ft3 = [ toProt(seq), toProt(seq[1:]), toProt(seq[2:]) ]

    curr_pos = 0
    best = []
    for p in framework:
        matches = [ (best_score(f, p, curr_pos), fn) for fn,f in enumerate(ft3) ]
        best.append(max(matches))
        curr_pos = best[-1][0][2]
    if best[0][0][0] < 90: continue
    if all(m[1] == 0 for m in best):
      print '+frame',
      f0 = ft3[0]
      fwk = []
      for (score, length, pos), frame in best:
        fwk.append(f0[pos:pos+length])
      cdr1 = f0[best[0][0][2]+best[0][0][1]:best[1][0][2]]
      cdr2 = f0[best[1][0][2]+best[1][0][1]:best[2][0][2]]
      cdr3 = f0[best[2][0][2]+best[2][0][1]:best[3][0][2]]
      print 'fwk', fwk[0], fwk[1], fwk[2], fwk[3],
      print 'cdr', cdr1, cdr2, cdr3
      print '+seq', f0
      sys.stdout.flush()
    else:
      print '-frame',
      fwk = []
      for (score, length, pos), frame in best:
        fwk.append(ft3[frame][pos:pos+length])
      print 'fwk', fwk[0], fwk[1], fwk[2], fwk[3],
      print (best[1][1] - best[0][1] + 3) % 3, (best[2][1] - best[1][1] + 3) % 3, (best[3][1] - best[2][1] + 3) % 3
