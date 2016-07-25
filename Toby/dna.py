import string

# correct ordering of NCBI tables
def remapCodonTable(table, src_order = 'TCAG', tgt_order = 'ACGT'):
  new_table = [None] * 64
  x = [tgt_order.index(src_order[_]) for _ in range(4) ]
  i = 0
  for a in range(4):
    for b in range(4):
      for c in range(4):
        new_table[(x[a] << 4) | (x[b] << 2) | x[c]] = table[i]
        i = i + 1
  return ''.join(new_table)

class CodonTable(object):
  trans_tbl = {
    # NCBI codon tables.
    1:  remapCodonTable('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    2:  remapCodonTable('FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG'),
    3:  remapCodonTable('FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    4:  remapCodonTable('FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    5:  remapCodonTable('FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG'),
    6:  remapCodonTable('FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    9:  remapCodonTable('FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'),
    10: remapCodonTable('FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    11: remapCodonTable('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    12: remapCodonTable('FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    13: remapCodonTable('FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG'),
    14: remapCodonTable('FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG'),
    15: remapCodonTable('FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    26: remapCodonTable('FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    21: remapCodonTable('FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG'),
    22: remapCodonTable('FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'),
    23: remapCodonTable('FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
  }

  def __init__(self, trans_tbl = 1):
    self.tbl = self.trans_tbl[trans_tbl]
    self.stop_codons = []
    for i in range(64):
      if self.tbl[i] == '*':
        self.stop_codons.append(self.uncodon(i))
    self.stop_codons = set(self.stop_codons)
    self.map = dict([ (self.uncodon(x), self.tbl[x]) for x in range(64) ])

  def codon(self, x):
    if type(x) in (type(0), type(0L)): return x
    if len(x) != 3: return None
    y = 0
    for _ in x:
      n = 'acgt'.find(_.lower())
      if n == -1: return None
      y = y * 4
      y = y + n
    return y

  def uncodon(self, x):
    return ''.join(['ACGT'[_&3] for _ in (x >> 4, x >> 2, x)])

  def translate(self, x):
    def c(codon, m = self.map):
      try:
        return m[codon.upper()]
      except:
        return 'X'
    return ''.join([ c(x[i:i+3]) for i in xrange(0, len(x), 3) ])

  def decode(self, x):
    try:
      return self.tbl[self.codon(x)]
    except TypeError:
      return None

  def synonymous(c1, c2):
    c1 = self.decode(c1)
    c2 = self.decode(c2)
    return c1 is not None and c1 == c2

_codon = CodonTable(trans_tbl = 1)

_s = 'XACMGRSVTWYHKDBN'
_enc_tab = [ -1 ] * 256
_dec_tab = [ 'X' ] * 256

for i in range(len(_s)):
  _enc_tab[ord(_s[i])] = i
  _enc_tab[ord(_s[i].lower())] = i
  _dec_tab[i] = _s[i]

def encode(x):
  return [ _enc_tab[ord(_)] for _ in x ]

def decode(x):
  return ''.join([ _dec_tab[_] for _ in x])

def revcomp(x, __trans = string.maketrans('XACMGRSVTWYHKDBNxacmgrsvtwyhkdbn', 'XTGKCYSBAWRDMHVNxtgkcysbawrdmhvn')):
  x = list(x)
  x.reverse()
  return string.translate(''.join(x), __trans)

def toProt(x, padded = 0):
  r = []
  for _ in range(0, len(x), 3):
    cod = x[_:_+3]
    COD = cod.upper()
    a = _codon.decode(COD)
    if a is not None and COD != cod:
      a = a.lower()
    if a is not None:
      r.append(a)
    else:
      r.append('X')
  if padded:
    return '.' + '..'.join(r) + '.'
  else:
    return ''.join(r)

# <codecell>

__all__ = [
  'revcomp',
  'toProt'
]
