#!env python
'''
This program loads a Fasta file and generates a frequency file that
looks just like those that are used by MEME and the MarkovBackground
class in EM.py

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
'''
import sys, re, os, math
from TAMO    import MotifTools  #Motif, top_nmers, oneletter, etc...
from TAMO.seq import Fasta



def main():
    seqsD = Fasta.load(sys.argv[1])
    seqs  = seqsD.values()
    for w in range(1,7):
        allnmers = permute(w)
        nmersT = MotifTools.top_nmers(w,seqs,'with counts','purge Ns')
        nmersD = {}
        total = 0
        for nmer in allnmers:
            nmersD[nmer] = 1 #Pseudo count
            total = total + 1
        for nmer,count in nmersT[:]:
            try: 
                rc = MotifTools.revcomplement(nmer)
                nmersD[nmer] = nmersD[nmer] + count
                nmersD[rc]   = nmersD[rc]   + count
                total = total + 2*count
            except KeyError:
                pass
        _t = nmersD.keys()
        _t.sort()
        print "# freq in %s (total %d with pseudocounts)"%(sys.argv[1],total)
        for nmer in _t:
            print "%-7s %20.17f"%(nmer,float(nmersD[nmer]) / total)
        sys.stdout.flush()

def permute(depth, letters=['A','C','G','T'], seqs=[''],curdepth=0):
    newseqs = []
    for seq in seqs:
        for letter in letters:
            newseqs.append(seq + letter)
    if depth > curdepth:
        return(permute(depth,letters,newseqs,curdepth + 1))
    else:
        return(seqs)



if __name__ == '__main__': main()
