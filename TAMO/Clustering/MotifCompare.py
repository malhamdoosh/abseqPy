#!env python
'''
This module contains low-level code for finding the optimal alignment between two motifs.
The central algorithm is called "minshortestoverhangdiff," which searches all alignments
of two motifs to find the one that minimizes the specificied distance/divergence metric.
The space of alignments is restricted by the minimum amount of overlap that is required
between motifs, which is in turn computed on the fly using the "OVLP" function.

Distance/divergence metrics include:
    -== To be completed ==-

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
'''
import sys, re, os, math, random

from  TAMO            import MotifTools
from  TAMO.util       import Arith
from  TAMO.MotifMetrics import ProbeSet

OVLP = lambda A,B,w=2: min(8,min(A.width,B.width)-w)
DFUNC = None
ACGT = list('ACGT')



def diffrange(self,other,Srange,Orange):
    '''Utility function: compute diff of '''
    '''self and other from self+Sstart, width of other'''
    POW     = math.pow
    Dtot    = 0
    divisor = float(len(Srange))
    
    '''Computes distance'''
    for si, oi in zip(Srange,Orange):
        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[si][L]) - POW(2,other.logP[oi][L]), 2 )
        col_dist  = math.sqrt(D)/math.sqrt(2.0)
        Dtot      = Dtot + col_dist
    _dist= Dtot/divisor
    #DBG
    #print self[min(Srange),max(Srange)+1], other[min(Orange),max(Orange)+1], _dist
    #print self.oneletter[min(Srange):max(Srange)+1]
    #print other.oneletter[min(Orange):max(Orange)+1], '%6.4f'%_dist
    #print 
    return _dist


def ENMrange(self,other,Srange,Orange):
    '''
    Utility function:  compute Expected Number of Mismatches
    '''
    '''Precompute Ps '''
    if (not  self.__dict__.has_key('P')) or (not  self.P):  self._computeP()
    if (not other.__dict__.has_key('P')) or (not other.P): other._computeP()

    '''Compute Overhangs'''
    sw = self.width
    ow = other.width
    shang = min(0,max(sw-Srange[-1]-1, sw-Srange[0]+1))
    ohang = min(0,max(ow-Orange[-1]-1, ow-Orange[0]+1))


    '''Compute Exp. # Mis'''
    Dtot = 0
    divisor = 1
    for si, oi in zip(Srange,Orange):
        S =  self.P[si]
        O = other.P[oi]
        D = 0
        for L in ACGT:
            D = D + S[L]*(1-O[L])
        Dtot      = Dtot + D
    _dist= Dtot/divisor
    return -_dist

def NLPO(A,B,thresh):
    p = probOvlp(A,B,thresh)
    return Arith.nlog10(p)

def NLBPO(A,B,thresh):
    p = probOvlpBinomial(A,B,thresh)
    return Arith.nlog10(p)

def probOvlpBinomial(A,B,thresh=0.7,verbose=None):
    if A.width >= B.width:
        Wide, Narrow = A, B
    else:
        Wide, Narrow = B, A

    RC = MotifTools.revcomplement
    newWide  = Wide[-1,Wide.width+1]
    if Wide.__dict__.has_key('bestWide'):
        bestWide = Wide.bestWide
    else:
        bestWideD = {}
        for x in newWide.bestseqs(thresh*newWide.maxscore):
            bestWideD[x[1]] = 1
        for x in bestWideD.keys():
            bestWideD[RC(x)] = 1
        Wide.bestWide = bestWideD.keys()
        bestWide = Wide.bestWide
    Wide = newWide

    D={}
    for i in range(len(bestWide)):
        D[i] = bestWide[i]
    P = ProbeSet(genome=D)
    matchNarrow = P.count_matching_probes(Narrow,thresh=thresh)
    
    if matchNarrow == 0:
        p = 1.0
        return p
    
    if not Narrow.__dict__.has_key('probNarrow'): Narrow.probNarrow = {}
    if Narrow.probNarrow.has_key(Wide.width):
        probNarrow = Narrow.probNarrow[Wide.width]
    else:
        probNarrow = estimate_frequency(Narrow,Wide.width,thresh=thresh)
        Narrow.probNarrow[Wide.width] = probNarrow

    p = Arith.binomialsumtail(probNarrow,len(bestWide),matchNarrow)
    print '\nD= %7.3f %9.4e %8d %7d %-14s %-20s %-14s %-20s'%(
        Arith.nlog10(p),probNarrow,len(bestWide),matchNarrow,
        A.family,A,B.family,B)


    return p


def probOvlp(A,B,thresh=0.7,verbose=None):
    if A.width >= B.width:
        Wide, Narrow = A, B
    else:
        Wide, Narrow = B, A

    RC = MotifTools.revcomplement
    if 1:
        newWide  = Wide[-1,Wide.width+1]
        if Wide.__dict__.has_key('bestWide'):
            bestWide = Wide.bestWide
        else:
            bestWideD = {}
            for x in newWide.bestseqs(thresh*newWide.maxscore):
                bestWideD[x] = 1
            for x in bestWideD.keys():
                bestWideD[RC(x)] = 1
            Wide.bestWide = bestWideD.keys()
            bestWide = Wide.bestWide
        Wide = newWide
    
        if Narrow.__dict__.has_key('bestNarrow'):
            bestNarrow = Narrow.bestNarrow
        else:
            bestNarrowD = {}
            for x in Narrow.bestseqs(thresh*Narrow.maxscore):
                bestNarrowD[x] = 1
            for x in bestNarrowD.keys():
                bestNarrowD[RC(x)] = 1
            bestNarrow = bestNarrowD.keys()
            Narrow.bestNarrow = bestNarrow
        
    #bestWide   = [x[1] for x in Wide.bestseqs  (thresh*Wide.maxscore)  ]
    #bestNarrow = [x[1] for x in Narrow.bestseqs(thresh*Narrow.maxscore)]

    countNarrow = len(bestNarrow)
    countWide   = len(bestWide)

    numtotal    = math.pow(4,Wide.width)
    fudgefactor = math.pow(4,Wide.width - Narrow.width)

    bestWideTups = [(x,MotifTools.revcomplement(x)) for x in bestWide]

    countBoth = 0
    for i in range(len(bestNarrow)):
        m_narrow = bestNarrow[i]
        delj = []

        for j in range(len(bestWideTups)):
            if (bestWideTups[j][0].find(m_narrow) >= 0) or (bestWideTups[j][1].find(m_narrow) >= 0):
                countBoth += 1
                delj.append(j)

        delj.reverse()  #Chew in from the back
        for j in delj:
            del(bestWideTups[j])


    if verbose: print '%10d %10d %10d %10d | %10d  %5d '%(
        countWide, numtotal, countNarrow *fudgefactor , countBoth , countNarrow, Wide.width - Narrow.width),
    
    p = Arith.hypgeomsummore(countWide,                 #Num Interesting
                             numtotal,                  #All k-mers
                             countNarrow * fudgefactor, #Number picked
                             countBoth                ) #Number found
    return p


def estimate_frequency(motif,k,samples=100000,thresh=0.7):
    #Build sequences
    estimate   = -30
    total      = 0
    totalcount = 0
    for i in range(40):
        long_string = 'ACGT'*(int(float(samples)*k/4))
        long_string = list(long_string)
        random.shuffle(long_string)
        random.shuffle(long_string)
        random.shuffle(long_string)
        long_string = ''.join(long_string)
        seqD = {}
        for i in range(samples):
            offset = k*i
            seqD[i] = long_string[offset:offset+k]
        P = ProbeSet(genome=seqD)
        count = P.count_matching_probes(motif,thresh=thresh)
        total      += float(samples)
        totalcount += float(count)
        f = totalcount/total
        d = math.fabs(f-estimate)/(estimate+0.00000001)
        estimate = f
        if d < 1e-4: break
        if i > 2 and totalcount > 100: break
        #print '%10d %10d %12.3e  %12.3e'%(totalcount, total, f, d)
    return estimate
    

def snegprobproductrange(self,other,Srange,Orange):
    '''sneg -- between 0 and 1?'''
    dS = -negprobproductrange(self,self,Srange,Srange)
    dO = -negprobproductrange(other,other,Orange,Orange)
    d  = negprobproductrange(self,other,Srange,Orange)
    return d /(max(dS,dO))
    
def negprobproductrange(self,other,Srange,Orange):
    '''Utility function: compute sum of products of probabilities

       SUM [
    '''
    LOG     = math.log
    '''Precompute Ps '''
    if (not  self.__dict__.has_key('P')) or (not  self.P):  self._computeP()
    if (not other.__dict__.has_key('P')) or (not other.P): other._computeP()

    '''Computes distance'''
    Dtot = 0
    divisor = 1
    for si, oi in zip(Srange,Orange):
        S =  self.P[si]
        O = other.P[oi]
        D = 0
        for L in ACGT:
            D = D + S[L]*O[L]
        Dtot      = Dtot + D
    _dist= Dtot/divisor
    return -_dist
    
def genNCBpenalty(penalty=0):
    return lambda s,o,Sr,Or,p=penalty: negcommonbitsrange(s,o,Sr,Or,PENALTY=p)

def NCBpenalty(penalty=0):
    d= negcommonbitsrange(self,other,Srange,Orange,PENALTY=0)

def snegcommonbitsrange(self,other,Srange,Orange):
    '''snegbits -- between 0 and 1?'''
    d= negcommonbitsrange(self,other,Srange,Orange,PENALTY=0)
    return 2+d/(max(self.width,other.width))

def negcommonbitsrange(self,other,Srange,Orange,PENALTY=0):
    '''Utility function: compute common bits in self and other

       Common bits are computed using this expression:
       SUM [ max(PENALTY,min(bitsA,bitsB)) ]            (Penalty usually == 0)

       where bitsA = Pb_ij * log(Pa_ij / Q_j)   [i = position, j = base letter]
    '''
    LOG     = math.log
    '''Precompute Ps '''
    if (not  self.__dict__.has_key('P')) or (not  self.P):  self._computeP()
    if (not other.__dict__.has_key('P')) or (not other.P): other._computeP()

    DIAGNOSTIC = 0

    commonbits = 0.0
    Q = self.background #Assume the same background for both motifs
    #for L in ACGT: Q[L] = 0.25
    #print Srange, Orange, self, other
    #print "O",other.width,len(other.P)
    #print "S", self.width,len( self.P)

    for si, oi in zip(Srange,Orange):
        if DIAGNOSTIC: print '%s %s'%(self.oneletter[si],other.oneletter[oi]),
        S =  self.P[si]
        O = other.P[oi]
        colbits = 0.0
        for L in ACGT:
            #Sbits    = S[L] * LOG(S[L]/Q[L])/LOG(2.)
            #Obits    = O[L] * LOG(O[L]/Q[L])/LOG(2.)
            #colbits += max(0,min(Sbits,Obits))
            Sbits    = O[L] * LOG(S[L]/Q[L])/LOG(2.)
            Obits    = S[L] * LOG(O[L]/Q[L])/LOG(2.)
            colbits += max(PENALTY,min(Sbits,Obits))
            if DIAGNOSTIC: print '%s %3.1f|%3.1f %5.2f|%5.2f  '%(
                L, S[L],O[L],Sbits,Obits),
        commonbits += max(0,colbits)
        if DIAGNOSTIC: print '%5.3f'%max(0,colbits), colbits

    if DIAGNOSTIC: print
    if DIAGNOSTIC: print '%', self.oneletter[min(Srange):max(Srange)+1]
    if DIAGNOSTIC: print ' ',other.oneletter[min(Orange):max(Orange)+1], '%4.2f'%commonbits
    if DIAGNOSTIC: print 

    return -commonbits

def entropyrange(self,other,Srange,Orange):
    '''Utility function: compute diff of '''
    '''self and other from self+Sstart, width of other'''
    LOG     = math.log
    Dtot    = 0

    '''Precompute Ps '''
    if (not  self.__dict__.has_key('P')) or (not  self.P):  self._computeP()
    if (not other.__dict__.has_key('P')) or (not other.P): other._computeP()
    '''Compute distance'''
    ds = []
    denom = 0
    for si, oi in zip(Srange,Orange):
        try:
            S =  self.P[si]
            O = other.P[oi]
        except:
            print self.oneletter,si,len(self.P),self.width
            print other.oneletter,oi,len(other.P),other.width
            sys.exit(1)
        D = 0
        for L in ACGT:
            s = S[L]
            o = O[L]
            D = D + (s*LOG(s/o) + o*LOG(o/s))/LOG(2.0)
            #print (s*LOG(s/o) + o*LOG(o/s))/LOG(2.0), 
        ds.append(D)
        #print
        Dtot      = Dtot + D
        try:
            denom = denom + self.bits[si] + other.bits[oi]
        except:
            print Srange
            print self.oneletter, si, len(self.bits), ['%3.1f'%x for x in self.bits]
            print Orange
            print other.oneletter, oi, len(other.bits), ['%3.1f'%x for x in other.bits]
            sys.exit(1)
         
    try: DIST = Dtot / denom
    except ZeroDivisionError: DIST = 1000
    print '%',self.oneletter[min(Srange):max(Srange)+1],
    print other.oneletter[min(Orange):max(Orange)+1], denom, '%6.4f'%DIST, '%6.4f'%Dtot, '%6.4f'%diffrange(self,other,Srange,Orange)
    return DIST


rcmemo = {}
def minshortestoverhangdiff(A,B,minoverlap=6,want_offset=None,DFUNC=None):
    if not DFUNC: DFUNC = diffrange
    if type(A) != type(B):
        print "Error: Attempted to compute alignment of objects that are not both Motifs"
        print "       types %s: %s  and %s: %s"%(M1,type(M1),M2,type(M2))
        sys.exit(1)
    if A.width < B.width:
        self  = B ;  other = A
    else:
        self  = A ;  other = B        
    wA = A.width
    wB = B.width
    wS = self.width
    wO = other.width
    if wS < wO: return 1

    lastS = wS - wO
    Dmin = 1000
    if want_offset: Ds = []
    key = '%s %s'%(self.oneletter,`[x['A'] for x in self.ll]`)
    if 1:  #MEMOize rc when doing many, many comparisons (factor x2.5 speedup)
        try: selfrc = rcmemo[key]
        except:
            selfrc = self.revcomp()
            rcmemo[key] = selfrc
    else: selfrc = self.revcomp()  #No Memo

    Ds = []
    if want_offset: offinfo = []
    for S,rc in [(self,0), (selfrc,1)]:
        for Ostart in range(max(0,wO-minoverlap)):
            Orange = range(Ostart,wO)
            Srange = range(0,len(Orange))
            D = DFUNC(S,other,Srange,Orange)
            if want_offset: offinfo.append((D,-Ostart,rc))
            Ds.append(D)
        Orange = range(0,wO)
        for Sstart in range(0,wS-wO+1):
            Srange = range(Sstart,Sstart+wO) 
            D = DFUNC(S,other,Srange,Orange)
            if want_offset: offinfo.append((D,Srange[0],rc))
            Ds.append(D)
        for ovlp   in range(wO-1,minoverlap,-1):
            Srange = range(wS-ovlp,wS)
            Orange = range(0,ovlp)
            D = DFUNC(S,other,Srange,Orange)
            if want_offset: offinfo.append((D,Srange[0],rc))
            Ds.append(D)
    if not want_offset:
        Ds.sort()
        return Ds[0]
    else:
        offinfo.sort()
        D,off,rc = offinfo[0]
        offset = off
        if rc and (self == A):  #We're not printing the RC of A, so adjust the offset
            offset = wA-wB-off
        elif      (self == B):  #We'll print the RC of B, so simply negate the offset
            offset = 0 - offset
        return offset, rc

def negcommonbitstest(t1,t2,OVLP_FCN=None):
    testdiff(t1,t2,OVLP_FCN,negcommonbitsrange)

def etestdiff(t1,t2,OVLP_FCN):
    testdiff(t1,t2,OVLP_FCN,entropyrange)

def testdiff(t1,t2,OVLP_FCN=None,DIFF_FCN=None):
    m1 = MotifTools.Motif_from_text(t1,bg=MotifTools.YEAST_BG)
    showdiffXvert(m1,t2,OVLP_FCN,DIFF_FCN)

def showdiffXvert(motif,seq,OVLP_FCN=None,DIFF_FCN=None):
    '''
    The funtion converts the sequence to a Motif, computes the D
    of the best alignment, and prints the alignment that generated
    that D.
    '''
    MSOdiff = minshortestoverhangdiff 
    if not OVLP_FCN:   OVLP_FCN = lambda A,B: min(min(A.width,B.width)-1,7)
    bg            = motif.background
    other         = MotifTools.Motif_from_text(seq,bg=bg)
    ovlp          = OVLP_FCN(motif,other)
    diff          = MSOdiff(motif,other,ovlp,DFUNC=DIFF_FCN)
    offset,rcflag = MSOdiff(motif,other,ovlp,'want_offset',DFUNC=DIFF_FCN)
    if rcflag: m=other.revcomp()
    else:  m=other
    print 'MSOdiff:  %8.4f %s%s%s'%(diff,' '*15,motif.oneletter,' '*(30-motif.width))
    print '          %8s %s%s%s'%(' ',' '*(15+offset),m.oneletter,' '*(30-offset-other.width))
    return diff

def averagetxts(txts,ovlp=2,VERBOSE=1):
    motifs = [MotifTools.Motif_from_text(x.strip()) for x in txts]
    return averagemotifs(motifs,ovlp,VERBOSE=VERBOSE)

def averagemotifs(motifs,ovlp=2,template=None,DFUNC=negcommonbitsrange,VERBOSE=1,prop=''):
    if not template: 
        Dmat = computeDmat(motifs)
        idx  = centroididx(Dmat)
        template = motifs[idx]

    for m in motifs:
        off, rc = minshortestoverhangdiff(template,m,OVLP(template,m),'want_offset',DFUNC=DFUNC)
        m.offset = off
        m.rc     = rc
        #Find most negative offset
    offsets = [m.offset for m in motifs]             ; offsets.sort()
    maxposs = [(m.offset + m.width) for m in motifs] ; maxposs.sort()
    minpos = -offsets[0]
    maxpos = maxposs[-1] + minpos
    pmotifs = []
    for m in motifs:
        if m.rc: _m = m.revcomp()
        else   : _m = m
        leftpad  = minpos + m.offset
        rightpad = maxpos - (leftpad + m.width)
        padded   = _m[-leftpad,_m.width+rightpad]
        #print '%s%s%s\t%s'%('*'*leftpad,_m.oneletter,'*'*rightpad,padded)
        pmotifs.append(padded)
    AVE = MotifTools.sum(pmotifs,[])
    if VERBOSE:
        for m in pmotifs:
            d = minshortestoverhangdiff(AVE,m,OVLP(AVE,m),DFUNC=DFUNC)
            print '%s   %5.3f'%(m.oneletter,d),
            if m.__dict__.has_key('key'): print m.key,
            if prop and m.__dict__.has_key(prop): print m.__dict__[prop],
            print
        print '-'*m.width
    return AVE


def coarsecluster(motifs):
    for i in range(len(motifs)):
        motifs[i].members = [motifs[i]]
        motifs[i].uniqid = i

    newmotifs = coalesce(motifs,-15) #Gather identical motifs
    nmotifs = len(newmotifs)
    nnewmotifs = 0
    while nnewmotifs != nmotifs:
        nmotifs    = len(newmotifs)
        newmotifs  = coalesce(newmotifs,-9)
        nnewmotifs = len(newmotifs)
        print 'ITERATION: ',nmotifs, nnewmotifs
    return newmotifs

def centroididx(Dmat):
    N = len(Dmat.keys())
    avedists = []
    for i in range(N):
        dtot = 0
        for j in range(N):
            if i == j: continue
            dtot += Dmat[i][j]
        avedists.append((dtot/N,i))
    avedists.sort()
    bestave,bestidx = avedists[0]
    return bestidx

def coalesce(motifs,DMIN=0.22):
    Dmat = computeDmat(motifs,'VERBOSE')
    nmotifs = len(motifs)
    claimed = []
    clusters = []
    for i in range(nmotifs):
        if i in claimed: continue
        cluster = [motifs[i]]
        for j in range(nmotifs):
            if i == j: continue
            if j in claimed: continue
            if Dmat[i][j] < DMIN:
                cluster.append(motifs[j])
        clusters.append(cluster)

    ans = []
    for _motifs in clusters:
        motifs = []
        for m in _motifs:
            motifs.extend(m.members)  #Add the members if motif is already a cluster
        miniDmat = computeDmat(motifs)
        centroid = motifs[centroididx(miniDmat)]
        AVE      = averagemotifs(motifs,ovlp=2,template=centroid)
        AVE.members = motifs
        ans.append(AVE)
    return ans
    
def computeDmat(motifs,VERBOSE=0,DFUNC=DFUNC):
    N = len(motifs)
    dmat = {}
    for i in range(N):                  #Initialize symmetric matrix
        dmat[i] = {}
        dmat[i][i] = 0.0
    if VERBOSE: print "              |%s|"%('-'*(N-1))  #Pretty status bar
    if VERBOSE: print "Computing ...  ",
    Nhalf = int(.293*N)
    for i in range(N-1):
        if VERBOSE:
            if i == Nhalf: sys.stdout.write('|'); sys.stdout.flush()
            else:          sys.stdout.write('.'); sys.stdout.flush()
        A = motifs[i]
        for j in range(i+1,N):
            B = motifs[j]
            #D = minshortestdiff(A,B)
            D = minshortestoverhangdiff(A,B,OVLP(A,B),DFUNC=DFUNC)
            dmat[i][j] = D              #Symmetric
            dmat[j][i] = D              #Symmetric
    if VERBOSE: print                               #End of Pretty status bar
    return dmat

def main():
    import pickle
    global DFUNC, DMIN
    DFUNC = negcommonbitsrange
    DMIN  = -9
    
    D = pickle.load(open(sys.argv[1]))
    for k,m in D.items(): m.key = k
    motifs = D.values()[0:40]
    results = coarsecluster(motifs)
    for m in results:
        print m.oneletter, len(m.members), ' '.join([x.key for x in m.members])

if __name__ == '__main__': main()
   
