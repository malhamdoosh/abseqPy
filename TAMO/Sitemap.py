#env python
#Copyright (2004) Whitehead Institute for Biomedical Research
#All Rights Reserved
#Author: David Benjamin Gordon

import sys, re, os, math, time, string, tempfile
from   TAMO     import MotifTools
from   TAMO.seq import Fasta

import getopt


def usage():
    print 'Usage: %s -f fasta_file -m motif_or_ambig '%(re.sub('^.*/','',sys.argv[0]))
    print '         [-t fraction (for PSSMs)]'
    print '         [-L <oneletterstring>]'
    print '         [-S <scale>   Currently 1/20 char/bases for pretty printing IGRs'
    print '\n'
    print 'Example: Find matches to ACCAT[CT] and CC...[AT][AT].GG, and'
    print '         label occurences of the first with "t" and the latter'
    print '         with "M"'
    print '\n         %s -f ypd_targets.fsa -m ACCATY,CCNNNWWNGG -L tM'%(sys.argv[0].split('/')[-1])
    print 

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:m:n:L:t:a:S:", ["help", "output="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    if not opts:
        usage()
        sys.exit(1)
        

    print "#" + ' '.join(sys.argv)
    fastafile, motiffile, motifnums, labels, thresh = (None, None, [], None, 0.7)
    ambigs = []

    scale   = 50.0 / 1000.0
    
    motifs = []
    for opt, value in opts:
        #print opt, value
        if   opt == '-f':  fastafile = value
        elif opt == '-m':  motifs.extend(MotifTools.txt2motifs(value))
        elif opt == '-n':  motifnums = [int(x) for x in value.split(',')]
        elif opt == '-L':  labels    = list(value)
        elif opt == '-t':  thresh    = float(value)
        elif opt == '-a':  ambigs.extend(value.split(','))
        elif opt == '-S':  scale     = float(value)
        
    probes = Fasta.load(fastafile)
    
    if motiffile:
        motifs.extend(TAMO.tamofile2motifs(motiffile))
    if ambigs:
        for ambig in ambigs:
            motifs.append( MotifTools.Motif_from_text(ambig,0.1) )
    if not motifnums:  motifnums = range(len(motifs))
    print '# %d: %s'%(len(motifs),motifnums)
    for i in range(len(motifnums)):
        motif = motifs[motifnums[i]]
        if labels and i < len(labels):
            txt = labels[i]
        else:
            txt = '%d'%i
        print '%-3s : %s %5.2f (%4.2f)'%(txt,motif,thresh*motif.maxscore,thresh)

    probehits = {}
    for key in probes.keys():
        hits_by_motif = []
        save_flag     = 0
        if re.search('[BDHU]',probes[key]): continue
        for num in motifnums:
            result = motifs[num].scan(probes[key],thresh*motif.maxscore)
            if result[0]:
                hits_by_motif.append(result)
                save_flag = 1
            else:
                hits_by_motif.append(None)
        if save_flag:
            probehits[key]=hits_by_motif

    #scale   = .1
    maxw = 40
    for key in probehits.keys():
        l       = len(probes[key])
        a       = list('-'* int(scale*l) )
        a.extend( list(' '*10 ) )
        desc    = []
        matches = probehits[key]
        for i in range(len(matches)):
            if matches[i]:
                subseqs,endpoints,scores = matches[i]
                for idx in range(len(subseqs)):
                    start,stop = endpoints[idx]
                    subseq     = subseqs[idx]
                    score      = scores[idx]
                    if labels and (i<len(labels)): ID = labels[i]
                    else                         : ID = '%d'%i
                    desc.append('%s %s %d-%d %4.2f '%(ID,subseq,start,stop,score))
                    start = int(start*scale)
                    for offset in range(10):
                        if a[start+offset] == '-':
                            if labels and (i < len(labels)):
                                a[start+offset] = labels[i]
                            else:
                                a[start+offset] = '%d'%i
                            break
        print '%-14s %s'%(key,''.join(a)),
        print ' '*max(0,maxw-len(a)), '| '.join(['%-27s'%x for x in desc])
        
    print
    print "Found matches in %d of %d input probes"%(len(probehits),len(probes))



def AmbigToRegExp(expression):
    expression=expression.replace('R','[RAG]')
    expression=expression.replace('Y','[YTC]')
    expression=expression.replace('W','[WTA]')
    expression=expression.replace('S','[SCG]')
    expression=expression.replace('M','[MAC]')
    expression=expression.replace('K','[KGT]')
    expression=expression.replace('H','[HATC]')
    expression=expression.replace('B','[BGCT]')
    expression=expression.replace('V','[VGAC]')
    expression=expression.replace('D','[DGAT]')
    expression=expression.replace('N','[NATCGRYWSMKHBVD-]')
    return expression




if __name__ == '__main__':
    main()
