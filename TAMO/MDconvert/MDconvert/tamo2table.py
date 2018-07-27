#!env python

#Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
#All Rights Reserved
#
#Author: David Benjamin Gordon

import sys, re, os, math, pickle, getopt
from TAMO              import MotifTools
from TAMO.util         import Arith
from TAMO.MotifCluster import minshortestdiff, rcmemo, minshortestoverhangdiff
from TAMO              import MotifMetrics
from TAMO.seq           import Fasta

GLOBALS = {}
GLOBALS['datafile']   = ''
GLOBALS['genomefile'] = 'YEAST'
GLOBALS['GIF']        = 1
GLOBALS['CONS']       = 0
GLOBALS['MAP']        = 1

def main():
    motifs = getarg('motifs')

    header()
    tablemotifs(motifs)
    footer()

def header():
    s =     ''
    s = s + '<html><head><title>MD Results</title></head>\n'
    s = s + '<body>\n'
    s = s + '<center><table border=1 style="text-align:center;"><tr bgcolor="lightblue">'
    s = s + '<td>Experiment</td>'
    s = s + '<td>#</td>'
    s = s + '<td>TRANSFAC matches to <br>to factor or to co-binding factors</td>'
    s = s + '<td>Motif</td>'
    s = s + '<td>#Targets</td>'
    s = s + '<td>%bound</td>'
    s = s + '<td>Enrichment</td>'
    s = s + '<td>ROC a.u.c.</td>'
    if getarg('MAP'):
        s = s + '<td>MAP/E</td>'
    if getarg('CONS'):
        s = s + '<td>%Conserved</td>'
        s = s + '<td>CI factor</td>'
    if getarg('PROGRAMS'):
        s = s + '<td>Programs</td>'
    s = s + '</tr>\n'
    print s
    sys.stdout.flush()

def footer():
    s =     ''
    s = s + '</table></center>\n'
    s = s + '</body></html>\n'
    print s
    sys.stdout.flush()

def tablemotifs(_motifs):
    file   = _motifs[0].file
    gifdir = 'logos_'+file+'.dir'
    
    motifs = [m.trimmed(0.25) for m in _motifs]
    cleanmotifs(motifs)

    top = getarg('top')
    if top: motifs = motifs[0:top]

    for m in motifs:
        if m.source:
            toks = m.source.split()
            m.name = toks[0].split('_')[0]
            if (len(toks) > 2) and ( toks[2] != '[]'):
                ids = toks[2].split(',')
                m.matchids = ', '.join([x.split('-')[-1] for x in ids])
            else: m.matchids ='&nbsp;'
        else:
            m.name = ''
            m.matchids = ''

    for i in range(len(motifs)):
        j = i + 1
        motifs[i].i = j
        if getarg('GIF'):
            gifname = motifs[i].giflogo(id='%s.%s'%(file,i), title=' ', scale=0.6, info_str=' ')
            if (not os.path.exists(gifdir)):   os.mkdir(gifdir)
            motifs[i].gifname = '%s/%s'%(gifdir,gifname)
            os.rename(gifname,motifs[i].gifname)

    txts = []
    for m in motifs:
        s=   ''
        s=s+ '<tr>'
        s=s+ '<td>%s</td>'%(m.name)
        s=s+ '<td>%d</td>'%(m.i) 
        s=s+ '<td>%s</td>'%(m.matchids)
        if getarg('GIF'): s=s+ '<td><img src="%s"></td>'%(m.gifname)
        else:   s=s+ '<td>%s</td>'%(m.oneletter)
        #s=s+ '<td>%5.2f</td>'%(m.i)
        try:
            s=s+ '<td>%3d</td><td>%5.1f%%</td>'%(m.numprobes,m.frac*100)
        except:
            s=s+ '<td>n/a</td><td>n/a</td>'
        s=s+ '<td>%5.2f</td>'%(Arith.nlog10(m.pvalue))
        try:
            s=s+ '<td>%5.2f</td>'%(m.ROC_auc)
        except:
            s=s+ '<td>%5.2f</td>'%(0.0)
        if getarg('MAP'):
            s=s+ '<td>%6.2f</td>'%(m.MAP)
        if getarg('CONS'):
            s=s+ '<td>%3.0f%% (%d)</td>'%(m.Cfrac*100,m.numconsmatchbound)
            s=s+ '<td>%d/%d -> %6.2f</td>'%(m.inC,m.inC+m.inU,m.CSimproved)
        if getarg('PROGRAMS'):
            s=s+ '<td>%s</td>'%("<br>".join(m.programs))
        s=s+ "</tr>\n"
        txts.append(s)
    print ''.join(txts);  sys.stdout.flush()

def cleanmotifs(motifs):
    for m in motifs:
        if not m.__dict__.has_key('programD'):
            m.programD = {}
        m.programs = m.programD.keys()
        m.programs.sort()


def usage(txt=''):
    if txt: print "Error: %s"%txt
    print 'Usage: %s -m motiffile '%(
        re.sub('^.*/','',sys.argv[0]))
    #print '          [-g|--genome genome.fsa]'
    #print '           [-w window (100)]'
    sys.exit(1)

def parse_opts():
    global GLOBALS
    short_opts = 'm:g:'
    long_opts  = ['genome=','top=']
    try:   opts, args = getopt.getopt(sys.argv[1:], short_opts, long_opts)
    except getopt.GetoptError:
        print getopt.GetoptError.__dict__
        usage()
    if not opts: usage()

    GLOBALS['args'] = args
    for opt,value in opts:
        if opt == '-m':                GLOBALS['motifs']     = MotifTools.txt2motifs(value)
        if opt in ['-g', '--genome']:  GLOBALS['genomefile'] = value
        if opt == '--top':             GLOBALS['top']        = int(value)

    #GLOBALS['genome'] = MotifMetrics.ProbeSet(getarg('genomefile'))
    #GLOBALS['data']   = Pvalues.Dataset(getarg('datafile'))

    

def getarg(varname):
    global GLOBALS
    if GLOBALS.has_key(varname):   return GLOBALS[varname]
    else:                          return None

if __name__ == '__main__':
    print "#" + ' '.join([x.replace(' ','\ ') for x in sys.argv])
    parse_opts()
    main()                                  

