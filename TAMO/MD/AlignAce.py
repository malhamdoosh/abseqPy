#!env python
"""
AlignAce.py -- Interface / Wrapper for the AlignACE program.

Download AlignACE separately from author's website (see TAMO/paths.py for instructions).
We have succesfully built AlignACE on linux and Win32 platforms (via cygwin)

This program is executable from the command line, and it runs the AlignACE several times
with different number seeds.  Type "$TAMO/MD/AlignAce.py" for options.

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
"""

import sys, re, os, math,  tempfile, random
import TAMO.paths
from   TAMO.MotifTools import Motif, print_motifs

def main():
    from TAMO import MotifMetrics
    if len(sys.argv) < 2:
        print "Usage: %s <fasta_file>"%(re.sub('^.*/','',sys.argv[0]))
        print '       [-w width (10)]    Model Width (note AlignACE allows gaps)'
        print '       [-iter    (10)]    Number of times to run AlignACE '
        print '       [-genome  fsafile] Genome (for computing background)'
        print '       [-gcback  (.38)    GC background (use 0.44 for human, 0.38 for yeast)]'
        sys.exit(1)

    print "#" + ' '.join([x.replace(' ','\ ') for x in sys.argv])

    fastafile = sys.argv[1]
    width     = 10
    valid_tfs = []
    iter      = 10
    genome    = 'YEAST'
    gcback    = 0.38
    
    for tok,i in zip(sys.argv,range(len(sys.argv))):
        if   tok == '-w'     : width = int(sys.argv[i+1])
        elif tok == '-valid' : valid_tfs.append(sys.argv[i+1])
        elif tok == '-iter'  : iter  = int(sys.argv[i+1])
        elif tok == '-gcback': gcback = float(sys.argv[i+1])
        elif tok == '-genome' :
            genome = sys.argv[i+1]
        elif tok == '-H250'  :
            genome = 'HUMAN_250'
            gcback = 0.44
        elif tok == '-Ch22'  :
            genome = 'Ch22'
            gcback = 0.44

    theMeta = MetaAce(fastafile,width,iter,gcback)

    Genome  = MotifMetrics.ProbeSet(genome)
    ids     = Genome.ids_from_file(fastafile)
    ids     = Genome.filter(ids)  #Only uses IDs that are actually in the Genome file

    motifs  = []
    motifs.extend(theMeta.results)

    for motif in motifs:
        motif.pvalue = Genome.p_value(motif,ids,factor=0.7)
        motif.church = Genome.church(motif,ids)
        for valid_tf in valid_tfs:
            motif.valid = Validate.validate(motif,valid_tf,'','Want Tuple')

    motifs.sort(lambda x,y: cmp(x.church,y.church))
    print_motifs(motifs,kmer_count=-1)
    

################################################################################
# AlignAce                                                                     #
# Execute, Load, and analyze AlignACE output                                   #
################################################################################
class AlignAce:
    '''
    ################################################################################
    # AlignAce                                                                     #
    # Execute, Load, and analyze AlignACE output                                   #
    ################################################################################

    Usage:
    theAce = AlignAce.AlignAce('filename.fsa',width=10)

    -or-

    theAce = AlignAce.AlginAce('filename.ace')

    The first form executes AlignACE with the specificied input sequences, and the second
    is used to parse alignace output from a previous run.

    Typically, the user is only interested in the ".motifs" member variable, which contains
    a list of all the discovered Motifs (as Motif objects).
    '''
    def __init__(self,file = '',width=10, args="",seed=''):
        self.EXE     = TAMO.paths.AlignACEdir + 'AlignACE'
        TAMO.paths.CHECK(self.EXE,'','AlignACE')
        self.lines   = []
        self.nmotifs = []
        self.motifs  = []
        self.outfile = ''
        self.width   = width
        self.AAargs  = args
        self.fastafile = ''
        self.outfile   = ''
        if seed: self.AAargs += ' -seed %d '%seed
        if (file[-4:] == '.fsa') or (file.find('.fasta')>0):  #Are we reading an AlignACE file or Analyzing a Fasta file?
            self.fastafile = file
        else:
            self.outfile = file
        if self.fastafile:
            self.outfile = re.sub('\.\w+$','',self.fastafile)
            self._execute(self.fastafile)
            self._parse()
        elif self.outfile:
            F = open(self.outfile,'r')
            self.lines=F.readlines()
            F.close()
            self._parse()
    def _execute(self,fastafile):
        command = '%s -i %s -numcols %d %s'%(self.EXE,fastafile,self.width,
                                             self.AAargs)
        FID = os.popen(command,'r')
        self.lines = FID.readlines()
        if FID.close():
            print 'Error executing command: \n\t\t%s'%command
    def _parse(self):
        i = 0
        while i < len(self.lines):
            line = self.lines[i]
            toks = line.split()
            if '-i' in toks and not self.fastafile:
                idx = toks.index('-i')
                self.fastafile = toks[idx+1]
            if (len(toks) > 0 and toks[0] == 'Motif'):
                seqs = []
                while 1:
                    i = i + 1
                    line = self.lines[i]
                    toks = line.split()
                    if toks[0][0] == '*': break
                    seqs.append(toks[0])
                M   = Motif(seqs)
                i = i + 1
                line = self.lines[i]
                toks = line.split()
                MAP = float(toks[2])
                if MAP > 1000: MAP = 0  #likely to be an AlignACE error
                M.MAP = MAP
                self.motifs.append(M)
            i = i + 1
        self.nmotifs = len(self.motifs)
    def __repr__(self):
        s = ''
        for M in self.motifs:
            s = s + "%s\t[%f]\n"%(M.__repr__(), M.MAP)
        return(s)

class MetaAce:
    '''
    Usage: theMetaAce = AlignAce.MetaAce(fastafile, width=10, iterations = 5, gcback=0.38)
    
    Class to invoke AlignACE many times and combine the results.  A random number
    generator is used to ensure that AlignACE is seeded differently each time.

    Typically, the user is only interested in the ".results" member variable, which contains
    a list of all the discovered Motifs (as Motif objects).
    '''
    def __init__(self,fastafile = '', width=10, iterations = 5, gcback=0.38):
        self.results = []
        self.seed    = int(random.random() * 1e9)
        self.gcback = gcback 
        self._execute(fastafile, width, iterations)
    def _execute(self,fastafile, width, iterations):
        v = 'verbose'
        for i in range(iterations):
            if v: print '#AlignACE iteration %d'%i,
            if v: sys.stdout.flush()
            A = AlignAce(fastafile,width,"-seed %d -gcback %f "%(self.seed,self.gcback))
            #print A.lines[6],
            self.results.extend(A.motifs)
            if v:
                print '# Iteration %2d: %3d motifs (%2d total)'%(
                    i,len(A.motifs),len(self.results))
                for m in A.motifs:
                    print "# %2d %s"%(i,m)
                sys.stdout.flush()
            self.seed = self.seed + 1
        if len(self.results) == 0: #Nothing found Retry recursively
            if (iterations > 1):
                self.seed = int(random.random() * 1e9)
                self._execute(fastafile, width, int(iterations/2)) 
            else:
                pass
                #os.system('cp %s tmpcur.fsa'%fastafile)
                #self.results.append(None)
        #print '>FOUND: %s<'%len(self.results),

    def flatten(self,thresh=0.5):
        '''
        flatten(thresh=0.5)
        This function does some messy reductions in metaACE lists, to cut the number of
        identical motifs reported.  A better solution is to use one of the clustering
        methods elsewhere in the TAMO package.
        '''
        motifs = self.results
        ans = []
        for m in motifs:
            match = 0
            for i in range(len(ans)):
                other = ans[i]
                if len(m) == len(other):
                    Df = m-other
                    Dr = m-other.revcomp()
                    if Df <= thresh  or Dr <= thresh: 
                        match = 1
                        if m.church == other.church:   #DBG 09-22-03 added to use enrichment
                            if m.pvalue < other.pvalue:
                                ans[i]=m
                        elif m.church < other.church:  #DBG 09-18-03 changed to < from >
                            ans[i]=m
                if match: break
            if not match: ans.append(m)
        ans.sort(lambda x,y: cmp(x.church,y.church)) #DBG 09-18-03 reversed sort order
        return ans
    
def flatten(motifs,thresh=0.5,prop='church',secondary_prop='pvalue'):
    '''
    flatten(motifs,thresh=0.5,prop='church',secondary_prop='pvalue'):
    This function does some messy reductions in metaACE lists, to cut the number of
    identical motifs reported.  A better solution is to use one of the clustering
    methods elsewhere in the TAMO package.

    "prop" and "secondary_prop" are used to settle ties between otherwise similar motifs.
    "church" refers to the AlignACE "group specificity score," and "pvalue" refers to
    the Enrichment score.
    '''
    ans = []
    for m in motifs:
        match = 0
        for i in range(len(ans)):
            other = ans[i]
            if len(m) == len(other):
                Df = m-other
                Dr = m-other.revcomp()
                if Df <= thresh  or Dr <= thresh: 
                    match = 1
                    if m.__dict__[prop] == other.__dict__[prop]:   #DBG 09-22-03 added to use enrichment
                        if m.__dict__[secondary_prop] < other.__dict__[secondary_prop]:
                            ans[i]=m
                    elif m.__dict__[prop] < other.__dict__[prop]:  #DBG 09-18-03 changed to < from >
                        ans[i]=m
            if match: break
        if not match: ans.append(m)
    ans.sort(lambda x,y,p1=prop,p2=secondary_prop: (cmp(x.__dict__[p1],y.__dict__[p1]) or  #DBG 09-18-03 reversed sort order
                                                    cmp(x.__dict__[p2],y.__dict__[p2])))
    return ans


if __name__ == '__main__': main()

