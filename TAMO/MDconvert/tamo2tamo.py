#!env python
#Copyright (2005) Whitehead Institute for Biomedical Research (except as noted below)
#All Rights Reserved
#
#Author: David Benjamin Gordon

import sys, re, os, math, time, string, tempfile, glob
from TAMO    import MotifTools
from TAMO.MD import AlignAce, Meme
from TAMO.util import Arith
from TAMO    import MotifMetrics
#from TAMO.altwebify import pick_genome

#PROBESET = MotifMetrics.ProbeSet('YEAST')
#PROBESET.factor = 0.7
probefile = None #Genome descriptor
PROBESET  = None #Genome object
letter    = 't'
fsafile   = None #Hand-entered fsa filep

def main():
    fsa_fcn = up_and_no_N

    parse()

    for filename in sys.argv[1:]:
        root = '.'.join(filename.split('.')[0:-1])
        ext  = filename.split('.')[-1]
        outname = '%s.%s%s'%(root,letter,ext)
        print '%-18s  --> %s'%(filename,outname)
        sys.stdout.flush()
        tamo2tamo(filename,outname)
        try: 
            pass
            #tamo2tamo(filename,outname)
        except:
            print "Error: Could not convert %s [[ %s ]]"%(
                filename, outname)
        


def parse():
    global probefile, PROBESET, letter, fsafile
    try:
        idx = sys.argv.index('-genome')
        del sys.argv[idx]
        probefile = sys.argv[idx]
        del sys.argv[idx]
        PROBESET = MotifMetrics.ProbeSet(probefile)
        PROBESET.factor = 0.7
    except: pass
    try:
        idx = sys.argv.index('-letter')
        del sys.argv[idx]
        letter = sys.argv[idx]
        del sys.argv[idx]
    except: pass
    try:
        idx = sys.argv.index('-f')
        del sys.argv[idx]
        fsafile = sys.argv[idx]
        del sys.argv[idx]
    except: pass
        
        

def tamo2tamo(file, outname):
    global probefile, PROBESET, fsafile
    
    motifs  = MotifTools.load(file)
    if fsafile:
        fsaname = fsafile
    else:
        fsaname = find_fsa(file)

    print '# FSA ',fsaname
    fsaD    = MotifMetrics.fasta2seqs(fsaname,'want_dict')
    probes  = fsaD.keys()
    if not probefile:
        PROBESET = MotifMetrics.ProbeSet('YEAST')
        #PROBESET= pick_genome(fsaname)
    #for key,seq in fsaD.items():
    #    PROBESET.probes[key] = seq

    print "# %d motifs"%len(motifs)
    for motif in motifs:
        #motif.pvalue, motif.church = 1,1  #Comment this!
        if motif.pvalue == 1: motif.pvalue = PROBESET.p_value(motif,probes,'v')
        if motif.church == 1: motif.church = PROBESET.church(motif,probes,'v')
        #if motif.E_site == None: motif.E_site = PROBESET.E_sitef(motif,probes,3,'v')
        #if motif.E_chi2 == None: motif.E_chi2 = PROBESET.E_chi2(motif,probes,None,'v')
        #if motif.E_seq  == None: motif.E_seq  = PROBESET.E_seq(motif,probes,'v')
        if motif.ROC_auc== None: motif.ROC_auc= PROBESET.ROC_AUC(motif,probes,'v')
        #if motif.MNCP   == None: motif.MNCP   = PROBESET.MNCP(motif,probes,'v')
        if motif.frac   == None: motif.frac   = PROBESET.frac(motif,probes,'v',0.7)
        if motif.numbound == 0:
            matching            = PROBESET.matching_ids(motif,[],factor=0.7)
            matchbound          = [x for x in matching if x in probes]
            motif.numbound      = len(probes)
            motif.nummotif      = len(matching)
            motif.numboundmotif = len(matchbound)
        if 0 and motif.CRA    == None:
            try:
                pass
                CRA, Cfrac = PROBESET.cons_ROC_AUC(motif,probes,'v',tuple='YES')
                motif.CRA = CRA
                motif.Cfrac = Cfrac
            except: pass
        
    MotifTools.save_motifs(motifs,outname)

def find_fsa(name,pathhint='../'):
    exists = os.path.exists
    root   = re.sub('\.\w*$','',name)
    smroot = re.sub('_.$'   ,'',root)
    tail   = root.split('/')[-1]
    parent = '/'.join(name.split('/')[:-2]) 

    if re.search('\.fsa$',name):
        if exists(name):
            return name
        elif exists(pathhint + name):
            return pathhint + name
    else:
        if exists(root + '.fsa'):
            return root + '.fsa'
        elif exists(parent + tail + '.fsa'):
            return parent + tail + '.fsa'
        elif exists(pathhint + root + '.fsa'):
            return pathhint + root + '.fsa'
        elif exists(smroot + '.fsa'):
            return smroot + '.fsa'
        elif exists(pathhint + smroot + '.fsa'):
            return pathhint + smroot + '.fsa'
    print '## ! Could not find fsa file for %s'%name
    return None
    

def up_and_no_N(name):  #Deprecated.... Use find_fsa instead
    root = re.sub('_N.ace','',name)
    ans = '../%s.fsa'%root
    return ans



if __name__ == '__main__': main()
