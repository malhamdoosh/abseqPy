#!env python
'''
This module uses a variant of the UPGMA algorithm for hierarchical clustering to
organize motifs into groups.  The "A" step (averaging) is performed by avergaing the
probability values in the aligned PSSMs within a cluster.

    1) Initialization
    1.1)  Assign each motif to a cluster
    1.2)  Compute Dmat of all clusters

    2) Iteration
    2.1)  Find the i and j with the smallest Distance
    2.2)  Create a new cluster (ij) which has n_i + n_j members
    2.3)  "Connect" i and j to (ij) and give each of the branchs D_ij/2 (better distance?)
    2.4)  Compute the distance from (ij) to all other clusters (except i and j)
    2.5)  Delete columns/rows in Dmat corresponding to i and j

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon

'''
import sys, re, os, math, getopt
GLOBALS = {}

from TAMO.Clustering.MotifCompare import *
from TAMO.Clustering              import MotifCompare
from TAMO                         import MotifTools

DFUNC = MotifCompare.DFUNC
DMAX  = 0.2

    
def main():
    motifs = getarg('motifs')
    print "Clustering %d motifs"%len(motifs)
    tree = UPGMA(motifs,DFUNC)
    print_tree(tree)
    motifs = slice_tree(tree,DMAX)
    for m in motifs:
        print m.oneletter,m.revcomp().oneletter,len(flatten_members(m))

def set_dfunc(_dfunc, _dmax):
    '''
    set_dfunc(dfunc, dmax)

    Set the distance/divergence/difference metric and threshold, and synchronize
    between this module and MotifCompare

    Examples:
    set_dfunc(MotifCompare.negcommonbitsrange, -8.0)
    set_dfunc(MotifCompare.diffrange, 0.23)
    '''

    global DFUNC, DMAX
    DFUNC = _dfunc
    MotifCompare.DFUNC = _dfunc
    DMAX  = _dmax

def UPGMAstep(motifs,DMAXrange,DFUNC=None):
    global DMAX
    tree = UPGMA(motifs,DFUNC)
    print_tree(tree)
    for dmax in DMAXrange:
        ans = slice_tree(tree,dmax)
        print '|--------- %4d clusters at %5.2f --------|'%(len(ans),dmax)
        for m in ans:
            print '%6d %-20s %-20s '%(
                len(flatten_members(m)),
                m.oneletter,m.revcomp().oneletter)
        print
        
def UPGMA(motifs,DFUNC=None):
    '''
    1) Initialization
    1.1)  Assign each motif to a cluster
    1.2)  Compute Dmat of all clusters

    2) Iteration
    2.1)  Find the i and j with the smallest Distance
    2.2)  Create a new cluster (ij) which has n_i + n_j members
    2.3)  "Connect" i and j to (ij) and give each of the branchs D_ij/2 (better distance?)
    2.4)  Compute the distance from (ij) to all other clusters (except i and j)
    2.5)  Delete columns/rows in Dmat corresponding to i and j
    '''
    #Now again.... with code...

    #1) Initialization
    #1.1)  Assign each motif to a cluster
    clusters = [c.copy() for c in motifs]
    for i in range(len(clusters)):
        clusters[i].members   = []
        clusters[i].idx       = i
        clusters[i].clustDave = 0.0
        clusters[i].clustDmax = 0.0
        clusters[i].clustDmin = 0.0

    #1.2)  Compute Dmat of all clusters
    Dmat = computeDmat(motifs,VERBOSE=1,DFUNC=DFUNC)
    SKIPD = {}

    #2) Iteration
    counter = len(clusters)-1

    print "            |%s|"%('-'*(len(clusters)-1))
    print "Clustering   ",
    while counter >= 1:
        sys.stdout.write('.')
        sys.stdout.flush()
        counter = counter - 1
        #2.1)  Find the i and j with the smallest Distance
        #2.5)  Delete columns/rows in Dmat corresponding to i and j
        i, j = find_min_ij(Dmat,SKIPD)
        SKIPD[i] = 1
        SKIPD[j] = 1
    
        #2.2)  Create a new cluster (ij) which has n_i + n_j members
        #2.3a) "Connect" i and j to (ij) 
        #2.3b) Give each of the branchs D_ij/2 (better distance?)
        newclust = collapse_clusters(clusters[i],clusters[j],Dmat,DFUNC)
        newclust.idx = len(clusters)
        clusters.append(newclust)
    
        #2.4)  Compute the distance from (ij) to all other clusters (except i and j)
        extendDmat(Dmat,clusters,SKIPD,DFUNC)

    for i in range(len(clusters)):
        c = clusters[i]
    tree = clusters[-1]
    #tree.dump = lambda x=tree:   print_tree(x)
    #tree.slice= lambda d,x=tree: slice_tree(x,d) 
    print 
    return tree

        

def find_min_ij(Dmat,SKIPD):
    idxs = [x for x in Dmat.keys() if not SKIPD.has_key(x)]
    N = len(Dmat)
    min = 1e+1000
    #print_dmat(Dmat,SKIPD)
    for i in idxs:
        for j in idxs:
            if i>=j: continue
            if Dmat[i][j] < min:
                imin, jmin, min = i,j,Dmat[i][j]
    return imin,jmin

def collapse_clusters(A,B,Dmat,_DFUNC=None):
    if not _DFUNC: _DFUNC = DFUNC
    A_members = flatten_members(A)
    B_members = flatten_members(B)
    motifs = A_members + B_members


    #Find centroid
    size = len(motifs)
    idxs = [x.idx for x in motifs]
    idxs.sort()
    avedists = []
    for i in idxs:
        dtot = 0
        for j in idxs:
            if i ==j: continue
            dtot += Dmat[i][j]
        avedists.append((dtot,i))
    avedists.sort()
    bestdist,bestidx = avedists[0]
    centroid = [x for x in motifs if (x.idx == bestidx)][0]
    #Monther's Hack
    #AVE = averagemotifs_counts
    AVE = averagemotifs(motifs,template = centroid,DFUNC=_DFUNC,VERBOSE=0)
    #Monther's Hack
    AVE.members = [A, B]

    dists = []
    dtot  = 0.0
    for m in motifs:
        D = minshortestoverhangdiff(m,AVE,OVLP(m,AVE),DFUNC=_DFUNC)
        dtot += D
        dists.append(D)

    AVE.clustDmin = min(dists)
    AVE.clustDmax = max(dists)
    AVE.clustDave = dtot / len(dists)
    
    return AVE

def flatten_members(M):
    if M.members:
        ans = []
        for m in M.members:
            ans.extend(flatten_members(m))
    else:
        ans = [M]
    return ans

def extendDmat(Dmat,clusters,SKIPD,_DFUNC=None):
    if _DFUNC==None: _DFUNC = DFUNC
    #print_dmat(Dmat)
    #Must assume that last cluster is the new one
    N = len(clusters)
    last = N-1
    lastM = clusters[last]
    Dmat[last] = {}
    for i in range(last):
        if SKIPD.has_key(i): continue
        A = clusters[i]
        D = minshortestoverhangdiff(A,lastM,OVLP(A,lastM),DFUNC=_DFUNC)
        Dmat[i][last] = D
        Dmat[last][i] = D
    #print_dmat(Dmat)


def print_tree(tree,level=0):
    if not tree: return
    txt = tree.oneletter  #This tree is a weakly-subclassed motif
    if tree.clustDmax > DMAX:  dtxt = ' '
    else:                      dtxt = '*'
    if not tree.members:  #This is a leaf
        leaf = '|-'
    else:
        leaf = '\_'
    ind = ' ' * (5 * level)
    info = ' ' * (5 * (20-level)-len(tree.oneletter)) + '%5.2f'%tree.clustDmax
    print ' %3d %s %s%s %s %s'%(level,ind,leaf,dtxt,txt,info)
    for subtree in tree.members:
        print_tree(subtree,level+1)

def create_tree_phylip(tree, level=0):
    if not tree: 
        return ''
    
    phylip = ''
    if not tree.members:
        phylip += tree.id + ':' + `tree.clustDmax` 
    else:
        phylip += '('
        for subtree in tree.members:
            phylip += create_tree_phylip(subtree, level+1) + ','    
        phylip = phylip.strip(',')
        if (level > 0):
            phylip += '):' + `tree.clustDmax`
        else:
            phylip += ');'
    return phylip
            
            
### added by Monther           
def print_tree_id(tree,level=0):
    if not tree: return
#     txt = tree.id  #This tree is a weakly-subclassed motif
    if tree.clustDmax > DMAX:  dtxt = ' '
    else:                      dtxt = '*'
    if not tree.members:  #This is a leaf
        leaf = '|-'
        txt = tree.id
    else:
        leaf = '\_'
        txt = 'Cluster (%5.6f)' % (tree.clustDmax)
    ind = '|  ' * level
    print ' %3d %s %s%s %s'%(level,ind,leaf,dtxt,txt)
    for subtree in tree.members:
        print_tree_id(subtree,level+1)


def slice_tree(tree,DMAX):
    ans = []
    if tree.clustDmax < DMAX: ans.append(tree) #This is a good rep of sub-branches
    elif not tree.members:    ans.append(tree) #This is a leaf
    else:
        for subtree in tree.members:
            subslice = slice_tree(subtree,DMAX)
            ans.extend(subslice)
    return ans

def print_dmat(D,SKIPD={}):
    keys = [x for x in D.keys() if not SKIPD.has_key(x)]
    keys.sort()
    for k in keys:
        print ' %4d ::'%(k),
        sk = [x for x in D[k].keys() if not SKIPD.has_key(x)]
        sk.sort()
        for j in sk:
            print '%4d: %6.2f  '%(j,D[k][j]),
        print
            
            
def usage(txt=''):
    '''
    Place information about command-line arguments here
    '''
    if txt: print "Error: %s"%txt
    print 'Usage: %s -m motifs'%(
        re.sub('^.*/','',sys.argv[0]))
    print ""

    print '         [-d threshold]        Distance threshold.  Default is 0.2'
    print '         [--dfunc <function>]  Distance metric.  Examples are NCB, and diffrange.  Any'
    print '                               "range" function in MotifCompare.py is acceptible'
    print ''
    print 'Some useful examples:'
    print ''
    print '   UPGMA.py -m <motiffile.tamo> -d  0.2  --dfunc diffrange (default)'
    print '   UPGMA.py -m <motiffile.tamo> -d -8.0 --dfunc NCB       (Count common bits, negated for minimization)'

    
    sys.exit(1)

def parse_opts():
    global GLOBALS
    global DFUNC, DMAX
    short_opts = 'm:'
    long_opts  = ['dfunc:']
    try:   opts, args = getopt.getopt(sys.argv[1:], short_opts, long_opts)
    except getopt.GetoptError:
        print getopt.GetoptError.__dict__
        usage()
    if not opts: usage()

    GLOBALS['args'] = args
    GLOBALS['motifs'] = []
    DFUNCtxt = None
    for opt,value in opts:
        if opt == '-m':                   GLOBALS['motifs'] = MotifTools.txt2motifs(value)
        if opt == '--dfunc':              DFUNCtxt = value
        if opt == '-d':                   DMAX     = float(value)

    # Deal with DFUNC and DMAX
    if DFUNCtxt == 'NCB':
        _DFUNC = MotifCompare.negcommonbits
    elif DFUNCtxt:
        try:
            exec ("_DFUNC = MotifCompare.%s"%DFUNCtxt)
        except:
            usage("No such distance metric: %s"%DFUNCtxt)
    if _DFUNC:  set_dfunc(_DFUNC,DMAX)

def getarg(varname):
    global GLOBALS
    if GLOBALS.has_key(varname):   return GLOBALS[varname]
    else:                          return None


if __name__ == '__main__': 
    parse_opts()
    print "#" + ' '.join([x.replace(' ','\ ') for x in sys.argv])
    main()                                  

