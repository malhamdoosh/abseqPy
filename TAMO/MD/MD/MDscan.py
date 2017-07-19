#!env python
"""
MDscan.py -- Interface / Wrapper for MDscan program.

Download MDscan separately from author's website (see TAMO/paths.py for instructions)
At this time, we believe that MDscan is only available on the linux platform

This program is executable from the command line, and it runs the MDscan several times
with different motif widths.  Type "$TAMO/MD/MDscan.py" for options.

Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
All Rights Reserved

Author: David Benjamin Gordon
"""
import sys, re, os, math, tempfile,  getopt
import TAMO.paths
from   TAMO            import MotifMetrics
from   TAMO.seq         import Fasta
from   TAMO.MotifTools import Motif,print_motifs

MDSCAN_DIR = TAMO.paths.MDscandir
MDSCAN_EXE = '%s/MDscan.linux'%MDSCAN_DIR


def usage(txt=''):
    if txt: print "Error: %s"%txt
    print 'Usage: %s -f fasta_file'%(
        re.sub('^.*/','',sys.argv[0]))
    print '           [--genome genomefile]'
    print '           [--range  start,stop]'
    print '           [--pcnt|--top  value]'
    print '           [--bgfile bgfile    ]'
    sys.exit(1)

def main():
    short_opts = 'f:'
    long_opts  = ['genome=', 'range=', 'top=', 'pcnt=', 'bgfile=']
    try:   opts, args = getopt.getopt(sys.argv[1:], short_opts, long_opts)
    except getopt.GetoptError:
        print getopt.GetoptError.__dict__
        usage()
    if not opts: usage()

    fastafile = ''
    top_count = 10
    top_pcnt  = None
    genome    = 'YEAST'
    w_start   = 8
    w_stop    = 15
    bgfile    = MDSCAN_DIR + 'yeast_int.bg'
    for opt,value in opts:
        if opt == '-f':         fastafile = value
        if opt == '--genome':   genome    = value
        if opt == '--top':      top_count = int(value)
        if opt == '--pcnt':     top_pcnt  = float(value)
        if opt == '--range':    w_start,w_stop= [int(x) for x in value.split(',')]

    print "#" + ' '.join(sys.argv)
    probeids = Fasta.keys(fastafile)
    Genome = MotifMetrics.ProbeSet(genome)

    probeids = Genome.filter(probeids)

    if top_pcnt: top_count = max(top_count,int(top_pcnt/100.0 * len(probeids)))

    theMeta = metaMDscan(fastafile,w_start,w_stop,top_count)

    for m in theMeta.motifs:
        m.pvalue = Genome.p_value(m,probeids,'v')
        m.church = Genome.church(m,probeids,'v')
        sys.stdout.flush()

    theMeta.motifs.sort(lambda x,y: cmp(x.pvalue,y.pvalue))
    print_motifs(theMeta.motifs)


class metaMDscan:
    def __init__(self,file='',width_start=6,width_stop=15,top_count=10,extra_args='',bgfile=''):
        self.results    = []
        self.motifs     = []
        self.file       = file
        self.width_start= width_start
        self.width_stop = width_stop
        self.top_count  = top_count
        self.extra_args = extra_args
        self.bgfile     = bgfile
        if file: self._execute()
        
    def _execute(self):
        TAMO.paths.CHECK(MDSCAN_EXE,'','MDscan')
        for width in range(self.width_start,self.width_stop):
            print "#MDscan Width iteration %3d (from %3d to %3d)"%(
                width,self.width_start,self.width_stop)
            sys.stdout.flush()
            result = MDscan(self.file, width,
                            self.extra_args + ' -t %d'%self.top_count,bgfile=self.bgfile)
            self.results.append(result)
            self.motifs.extend (result.motifs)
            for m in result.motifs:
                print '# ',m.summary()
        

 
class MDscan:
    '''
    Class for encapsulating (and processing) MDscan and its output
    '''
    def __init__(self,file='',width=10,extra_args='',bgfile=''):
        self.fsafile = ''
        self.outfile = ''
        self.motifs  = []
        self.bgfile  = bgfile
        if (file[-4:] == '.fsa') or (file.find('.fasta')>0):
            self.fsafile = file
        else:
            self.outfile = file
        self.extra_args = extra_args
        self.width      = width

        if self.fsafile:
            self._execute()
            self._parse()
        elif self.outfile:
            F = open(self.outfile,'r')
            self.lines=F.readlines()
            F.close()
            self._parse()

    def _execute(self):
        TAMO.paths.CHECK(MDSCAN_EXE,'','MDscan')
        if not self.bgfile: bgfile = MDSCAN_DIR + 'yeast_int.bg'
        else:               bgfile = self.bgfile
        command = '%s -i %s -f %s -w %d -s 30 -r 5 %s -g 1'%(
            MDSCAN_EXE, self.fsafile, bgfile,
            self.width, self.extra_args)
        print '#',command
        FID = os.popen('{ %s ; } 2>&1'%command,'r')
        self.lines = [x.strip() for x in FID.readlines()]
        if FID.close():
            print 'Error executing command: \n\t\t%s'%(command)
            for line in self.lines: print line
    def _parse(self):
        'Parse MDscan file'
        alloutput = '\n'.join(self.lines)
        premotifs = alloutput.split('\nMtf ')
        print len(premotifs)
        for pm in premotifs:
            sublines = pm.split('\n')
            score, seednum = 0,0
            seqs = []
            for line in sublines:
                if line.find('Final Motif') == 0:
                    toks    = line.split()
                    score   = float(toks[6])
                    seednum = int(toks[8])
                if line.find('>') == 0:
                    seqs.append(line.split()[-1])
            #print "SEQS: ",seqs
            if seqs:
                m = Motif(seqs)
                m.MAP = score
                m.seednum = seednum
                self.motifs.append(m)
        

def genomebg(infile,outfile):
    EXE = MDSCAN_DIR + 'genomebg.linux'
    fsaD   = Fasta.load(infile)
    tmpfsa = tempfile.mktemp()
    Fasta.write(fsaD,tmpfsa,linelen=1000000000)
    CMD = '%s -i %s -o %s'%(EXE,tmpfsa,outfile)
    FID = os.popen('( %s ;) 2>&1'%CMD,'r')
    for line in FID.readlines(): print line
    if FID.close(): print "Exited"
    os.unlink(tmpfsa)
                                  
if __name__ == '__main__': main()                                  
    
