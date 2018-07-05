#!env python
"""
MDscan.py -- Interface / Wrapper for MDscan program.

Download MEME separately from author's website (see TAMO/paths.py for instructions)

This program is executable from the command line, and it runs the MEME with a
set of default parameters.   Type "$TAMO/MD/Meme.py" for options.

(There aren't too many options at this point.  This module is most useful for reading
 output from previous MEME runs.)

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
"""
import sys, re, os
import TAMO.paths
from   TAMO            import MotifMetrics
from   TAMO.MotifTools import * 

def main():
    if len(sys.argv) < 2:
        print "Usage: %s <fasta_file>"%(re.sub('^.*/','',sys.argv[0]))
        print "     [-genome genomefile.fsa]   Genome file (for computing Enrichment, etc..."
        print "     [-bfile  file          ]   File for Markov Background Model"
        print '     [-bigdata              ]   Adds "-maxsize 2000000" for large datasets'
        sys.exit(1)

    fastafile = sys.argv[1]
    width     = 0
    valid_tfs = []
    iter      = 10
    genome    = 'YEAST'
    xtra      = ''
    bfile     = None
    
    for tok,i in zip(sys.argv,range(len(sys.argv))):
        if   tok == '-w'     : width = int(sys.argv[i+1])
        elif tok == '-human' : genome = 'HUMAN'
        elif tok == '-H250'  : genome = 'HUMAN_250'
        elif tok == '-Ch22'  : genome = 'Ch22'
        elif tok == '-genome': genome = sys.argv[i+1]
        elif tok == '-bigdata': xtra  = '-maxsize 2000000'
        elif tok == '-bfile':  bfile = sys.argv[i+1]

    theMeme = Meme(fastafile,width,xtra,genome,bfile)
    Genome  = MotifMetrics.ProbeSet(genome)
    ids     = theMeme.probes
    #ids     = Genome.ids_from_file(fastafile)
    
    motifs = theMeme.motifs
    for motif in motifs:
        motif.pvalue = Genome.p_value(motif,ids,'v')
        for valid_tf in valid_tfs:
            motif.valid = Validate.validate(motif,valid_tf,'','Want Tuple')

    print_motifs(motifs)

    print '#'*80
    for line in theMeme.lines:
        print line,

class Meme:
    '''
    Class for encapsulating (and processing) MEME output
    '''
    def __init__(self,file='', width='', extra_args='',genome='YEAST',bfile=None):
        self.EXE    = TAMO.paths.MEMEdir + 'bin/meme'
        TAMO.paths.CHECK(self.EXE,'','MEME')
        if not bfile:
            if genome == 'YEAST':
                bfile = TAMO.paths.MEMEdir + 'tests/yeast.nc.6.freq'
            else:
                print "Need to specify a background frequency file if not yeast."
                print "Perhaps something like -bfile %s"%(
                    TAMO.paths.MEMEdir + 'tests/yeast.nc.6.freq')
                sys.exit(1)
        self.extra_args = extra_args
        self.lines  = []
        self.motifs = []
        self.probes = []
        self.fastafile = ''
        self.outfile   = ''
        self.args   = '-dna -revcomp -nmotifs 1 -mod zoops -nostatus -text -maxsize 200000 ' 
        if bfile: self.args = self.args + '-bfile %s '%bfile 
        self.width  = '-minw 6 -maxw 18 '
        if width:
            self.width = '-w %d '%width
        if (file[-4:] == '.fsa') or (file.find('.fasta')>0):
            self.fastafile = file
        else:
            self.outfile = file
        if self.fastafile:
            self._execute()
            self._parse()
        elif self.outfile:
            F = open(self.outfile,'r')
            self.lines=F.readlines()
            F.close()
            self._parse()
    def _execute(self):
        command = '%s %s %s %s %s'%(self.EXE,    self.fastafile,
                                    self.width,  self.args,
                                    self.extra_args)
        print '#',command
        FID = os.popen(command,'r')
        #self.lines = FID.readlines()
        while 1:
            line = FID.readline()
            if line == '': break
            self.lines.append(line)
            sys.stdout.write(line)
            sys.stdout.flush()
        if FID.close():
            print 'Error executing command: \n\t\t%s'%(command)
            for line in self.lines: print line,
    def _parse(self):
        'Parse MEME file'
        i = 0
        num = 0
        background = {}
        for i in range(len(self.lines)):
            toks = self.lines[i].split()
            for tok in toks: tok.strip()
            if len(toks) == 2 and not self.fastafile and toks[0] == 'DATAFILE=':
                self.fastafile = toks[1]
            if len(toks) < 4: continue
            if toks[0] == 'Background' and toks[2] == 'frequencies':
                toks = self.lines[i+1].split()
                for j in [0,2,4,6]:
                    background[toks[j]] = float(toks[j+1])
            #print self.lines[i],
            if toks[0] == "MOTIF":
                m = Motif([],background)
                num    = int(toks[1])
                m.evalue = float(toks[13])
            elif (toks[0] == 'Motif' and int(toks[1]) == num and
                  toks[3] == 'BLOCKS'):
                offset = 3
                j      = 0
                while self.lines[i+offset+j][0:2] != '//':
                    #print self.lines[i+offset+j],
                    _toks = self.lines[i+offset+j].split()
                    seq = _toks[-2]
                    seq = re.sub('X','N',seq)
                    m.seqs.append(seq)
                    j = j + 1                
            elif (toks[0] == 'Motif' and int(toks[1]) == num and
                  toks[3] == 'probability'):
                offset = 3
                j      = 0
                counts = []
                while self.lines[i+offset+j][0:10] != '------------'[0:10]:
                    probs = self.lines[i+offset+j].split()
                    _d = {}
                    for key,val in zip(['A', 'C', 'G', 'T'], probs):
                        _d[key] = float(val)
                    counts.append(_d)
                    j = j + 1
                for c in counts:
                    print c
                m.compute_from_counts(counts,0.00001)
                self.motifs.append(m)
            elif (toks[0] == 'Sequence' and toks[1] == 'name' and
                  toks[2] == 'Weight'):
                offset = 2
                j      = 0
                while self.lines[i+offset+j][0:5] != '*****':
                    _toks = self.lines[i+offset+j].split()
                    self.probes.append(_toks[0])
                    if len(_toks) > 5:
                        self.probes.append(_toks[3])
                    j = j + 1

if __name__ == '__main__': main()
