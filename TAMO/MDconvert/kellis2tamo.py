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
probefile = None
PROBESET  = None

def main():
    fsa_fcn = up_and_no_N

    parse()

    FID = open(sys.argv[1])
    tokss = [x.strip().split(',') for x in FID.readlines()]
    FID.close()

    D = {}
    for expt,motif,score,source in tokss:
        print expt,motif
        if expt == 'Category': continue
        if motif == 'x': continue
        motif = MotifTools.Motif_from_text(motif)
        motif.kellis = float(score)
        motif.source = source
        try: D[expt].append(motif)
        except: D[expt] = [motif]

    for expt,motifs in D.items():
        root = expt
        ext  = 'cons'
        if root[0:3] == 'Rnd':
            num = re.sub('.*_','',root)
            if len(num) == 1:
                root = re.sub('_','_00',root)
            else:
                root = re.sub('_','_0',root)
            root = re.sub('Rnd','random_',root)
        outname = '%s.t%s'%(root,ext)
        print '%-18s  --> %s'%(root,outname)
        sys.stdout.flush()
        motifs2tamo(motifs,outname)
        try: 
            pass
            #tamo2tamo(filename,outname)
        except:
            print "Error: Could not convert %s [[ %s ]]"%(
                filename, outname)
        


def parse():
    global probefile, PROBESET
    try:
        idx = sys.argv.index('-genome')
        del sys.argv[idx]
        probefile = sys.argv[idx]
        del sys.argv[idx]
        PROBESET = MotifMetrics.ProbeSet(probefile)
        PROBESET.factor = 0.7
    except: pass
        

def motifs2tamo(motifs, outname):
    global probefile, PROBESET
    
    fsaname = find_fsa(outname)
    fsaD    = MotifMetrics.fasta2seqs(fsaname,'want_dict')
    probes  = fsaD.keys()
    if not probefile:
        PROBESET = MotifMetrics.ProbeSet('YEAST')
        #PROBESET= pick_genome(fsaname)
    #for key,seq in fsaD.items():
    #    PROBESET.probes[key] = seq

    print "# %d motifs"%len(motifs)
    for motif in motifs:
        if motif.pvalue == 1: motif.pvalue = PROBESET.p_value(motif,probes,'v')
        if motif.church == 1: motif.church = PROBESET.church(motif,probes,'v')
        if motif.E_site == None: motif.E_site = PROBESET.E_sitef(motif,probes,3,'v')
        #if motif.E_chi2 == None: motif.E_chi2 = PROBESET.E_chi2(motif,probes,None,'v')
        if motif.E_seq  == None: motif.E_seq  = PROBESET.E_seq(motif,probes,'v')
        if motif.ROC_auc== None: motif.ROC_auc= PROBESET.ROC_AUC(motif,probes,'v')
        if motif.MNCP   == None: motif.MNCP   = PROBESET.MNCP(motif,probes,'v')
    MotifTools.save_motifs(motifs,outname)

def find_fsa(name,pathhint='../'):
    exists = os.path.exists
    root   = re.sub('\.\w*$','',name)
    smroot = re.sub('_.$'   ,'',root)

    print root
    if re.search('\.fsa$',name):
        if exists(name):
            return name
        elif exists(pathhint + name):
            return pathhint + name
    else:
        if exists(root + '.fsa'):
            return root + '.fsa'
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
