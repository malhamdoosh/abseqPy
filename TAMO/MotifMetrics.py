#!/usr/bin/env python
'''
MotifMetrics -- Kitchen sink of routines for evaluating motifs on genomes

This module is used both as a command-line script and as a repository
for motif metrics.  For a summary of the command-line usage, just type

    python MotifMetrics.py

for help.

-= Summary =-

The MotifMetrics module performs fast motif evaluation by converting
python representations of motifs (PSSMs or regular expressions) and
DNA sequences into efficient C++ data structures.

-= Motif Metrics =-

Motif metrics represent different ways to assign a score to the number
of binding sites that match a motif model that can be found within a
subset of genomic sequences (such as those identified by microarray or
functional annotation) as compared to the occurences in the entire
genome (or at least the entire sequence represented in the microrray.)
Historically, these metrics have generally been limited to intergenic
sequences, but they can be extended to complete sequences as
microarray technology progresses.

This module provides an interface to more than a dozen different
metrics.  These include:

   The "Enrichment" metric used in Harbison et al. Nature 2004
   The "Group Specificity Metric" used by AlignACE (Hughes et al.)
   The "ROC auc" and "MNCP" metrics described by Clarke et al.
   and several other variants and experimental metrics.  

Although many metrics are based on approaches from probability theory,
the scores themselves may not always be interpreted directly in terms
of significance, particularly if 1) the motif is learned from the same
data on from which the score(s) is computed or 2) the motif is one of
many motifs being evaluated on the same data.

-= Input =-

The metrics used in the module (generally) require these items as
input:

   1) A motif (specified as a regular expression, a IUPAC ambiguity
   code, or as a "tamo"-formatted PSSM).
   
   2) A genome (Usually specified as a fasta-file with individual
   entries representing each intergenic region or UTR region.)

   3) A list (or fasta-formatted file) of "interesting" genes or
   intergenic regions, using the same names as in the genome
   fasta-formatted file.

For example, to compute the "Enrichment" metric of the GCN4 motif
among intergenic probes identified in a GCN4 ChIP-chip assay in yeast,
the user must supply:

   1) The motif ("TGASTCA")

   2) A genome (say, the probe sequences on the intergenic array from
   the ChIP experiment: "yeast6k_probes.fsa")

   3) A list of the "bound" probes.  (Say, "GCN4_YPD.fsa")

NOTE: The first time a fasta-formated file is loaded into a
MotifMetrics.ProbeSet object (often given the variable name "genome"),
the object is cached to disk.  This means that the first user to load
a fasta-formatted file needs to have write permission in the same
directory.

-= Quick and Dirty motif discvery =-

Given a motif metric, it is possible to evaluate many candidate motifs
to find the ones with the best scores.  This is essentially a naive
motif discvery, that is limited only by the available computer time.
One built-in form of this available with the "-n" option, which
computes the specified metric for all k-mers (without ambiguity codes)
in the set of "bound" or "differentially regulated" sequences.  On a
P4-2000, about a thousand such motifs can be evaluated in about five
numutes.
   

-= Examples =-

At the command line:

% MotifMetrics.py GCN4_YPD.fsa -g $TAMOdata/yeast/yesat6k_probes.fsa TGASGCA

In a python file:

Genome = MotifMetrics.ProbeSet('../yeast/yeast6k_probes.fsa')
ids    = Fasta.ids('GCN4_YPD.fsa')
e      = Genome.Enrichment('TGASTCA',ids)

print e, -math.log(e)/math.log(10)

Copyright (2005) Whitehead Institute for Biomedical Research
All Rights Reserved
Author: David Benjamin Gordon

'''
#Python imports
import sys,os,re
import string
import shelve
import pickle
import math
#from Stats    import stats   #If you've got Stats, uncomment for KS and chi-square tests

import TAMO.paths

from TAMO.seq  import Fasta
from TAMO.util import Arith
from TAMO      import MotifTools
from TAMO.MD   import MDsupport
from TAMO.util.PermuteTools import gapped_words
from TAMO.util.PermuteTools import dimer_groups as dimer_words

def main():
    genome   = 'YEAST'
    motifs   = []
    N_file   = ''
    N_genome = ''
    toMotif  = ''
    tamofile = ''
    tamoidx  = -1
    tamoidxs = []
    factor   = 0.9
    find_best= 0
    gapped   = 0 
    TRANSFAC = 0
    TFMAT    = None  #Not Used
    verbose  = 0
    reducefile=None
    subfsa   = None
    bysitek  = 2
    byfactk  = 2
    omits    = []
    
    method   = 'E' #Could be "spec", "bysite", or "chi2"
    #spec     = None
    #bysite   = None
    #chi2     = None

    if len(sys.argv) < 2:
        print 'Usage: %s [idfile | fastafile] <options> <ambiguity code> <ambiguity code> ...'%(re.sub('^.*/','',sys.argv[0]))
        print '                --- Default is Yeast Probe Set ---'
        print '   -bysite         Count sites instead of probes   (motifs only, or with -M)'
        print '   -chi2           Chi2 over distribution of sites (motifs only, or with -M)'
        print '   -spec           Specificity instead of Enrichment'
        print '   -binomial       Binomial Enrichmet'
        print '' 
        print '   -y              Yeast Genome <default>'
#       print '   -Y              Yeast Consensus Genome'
#       print '   -Y2K, -Y5C      Yeast Upstream 2000 bp, 500 bp'
#       print '   -H              Human Probe set'
#       print '   -H 250 -H 1000  Expanded Human Probe set'
#       print '   -B              Bacterial Orf SLR14.2'
#       print '   -Ch22           Human Chr 22, according to Martone et. al'
        print '   -g|-genome file Genome specified as fsa file'
        print ''
        print '   -P              Permute -- all n-mers'
        print '   -n width        use top n-mers in file'
        print '   -N width        use top n-mers in genome database'
        print '   -M              generate PSSMs from text'
        print '   -t file <#>     get motif from TAMO-format file'
        print '   -m file:#       get motif from TAMO-format file or from pickle of motif list or dict'
        print ''
        print '   -V              Verbose'
        print '   -f frac         cutoff factor to multiply max score 0.9 default'
        print '   -X              for PSSMs, find best cutoff'
#       print '   -TF             Transfac PSSMs'  #Transfac not included in this distribution
        print '   -O filename     Omit IDs in file'
        print '   -R filename     Reduce genome to probes in file'
        print '   -sub            Replace "genome" entries with those in fsafile'
        print ''
        
        print '   -gap            test gapped motifs, gap 0-12'
        print '   -syl size       Look for sylA-gap-sylB words'
        print '   -dimer size     Look for sylA-gap-sylA (and revcomp) words'
        sys.exit(1)

    print "#" + ' '.join(map(lambda x: re.sub(' ','\ ',x), sys.argv))

    if '-y' in sys.argv:
        idx = sys.argv.index('-y')
        del sys.argv[idx]
        genome = 'YEAST'

    if '-Y' in sys.argv:
        idx = sys.argv.index('-Y')
        del sys.argv[idx]
        genome = 'YEAST_CONSERVED'

    if '-Y2K' in sys.argv:
        idx = sys.argv.index('-Y2K')
        del sys.argv[idx]
        genome = 'YEAST_2000_UP'

    if '-Y5C' in sys.argv:
        idx = sys.argv.index('-Y5C')
        del sys.argv[idx]
        genome = 'YEAST_500_UP'

    if '-B' in sys.argv:
        idx = sys.argv.index('-B')
        del sys.argv[idx]
        genome = 'BAC_ORF'

    if '-Ch22' in sys.argv:
        idx = sys.argv.index('-Ch22')
        del sys.argv[idx]
        genome = 'Ch22'

    if '-genome' in sys.argv:
        idx = sys.argv.index('-genome')
        del sys.argv[idx]
        genome = sys.argv[idx]
        del sys.argv[idx]

    if '-g' in sys.argv:
        idx = sys.argv.index('-g')
        del sys.argv[idx]
        genome = sys.argv[idx]
        del sys.argv[idx]

    if '-H' in sys.argv:
        idx = sys.argv.index('-H')
        if len(sys.argv) < (idx + 2) or  not sys.argv[idx+1].isdigit():
            genome = 'HUMAN'
        else:
            s   = sys.argv[idx+1]
            genome = 'HUMAN_%s'%s
            del sys.argv[idx]
        del sys.argv[idx]

    if '-P' in sys.argv:
        idx = sys.argv.index('-P')
        del sys.argv[idx]
        motifs.extend(permute(['A','C','T','G'], 7))

    if '-n' in sys.argv:
        idx    = sys.argv.index('-n')
        N_file = int(sys.argv[idx+1])
        del sys.argv[idx]
        del sys.argv[idx]
        print "#Using top %d-mers"%N_file

    if '-N' in sys.argv:
        idx      = sys.argv.index('-N')
        N_genome = int(sys.argv[idx+1])
        del sys.argv[idx]
        del sys.argv[idx]
        print "#Using top %d-mers"%N_genome

    if '-M' in sys.argv:
        idx      = sys.argv.index('-M')
        del sys.argv[idx]
        print "#Forming PSSMs from text"
        toMotif = 1

    if '-R' in sys.argv:
        idx    = sys.argv.index('-R')
        reducefile = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]
        print "#Reducing probeset according to %s"%reducefile

    if '-O' in sys.argv:
        idx    = sys.argv.index('-O')
        _t     = sys.argv[idx+1]
        omits  = Fasta.keys(_t)
        del sys.argv[idx]
        del sys.argv[idx]
        print "#Omiting probes from  %s"%_t

    if '-gap' in sys.argv:
        idx = sys.argv.index('-gap')
        gapped = 1
        del sys.argv[idx]
        motifs.extend(gapped_motifs(sys.argv[-1],24))

    if '-syl' in sys.argv:
        idx = sys.argv.index('-syl')
        gaplength = int(sys.argv[idx+1])
        motifs.extend(gapped_words(int(sys.argv[idx+1])))
        del sys.argv[idx]
        del sys.argv[idx]

    if '-dimer' in sys.argv:
        idx = sys.argv.index('-dimer')
        gaplength = int(sys.argv[idx+1])
        motifs.extend(dimer_words(int(sys.argv[idx+1])))
        del sys.argv[idx]
        del sys.argv[idx]

    if '-m' in sys.argv:
        idx     = sys.argv.index('-m')
        filetxt = sys.argv[idx+1].split(':')
        motifs.extend(MotifTools.txt2motifs(filetxt))
        del sys.argv[idx]
        del sys.argv[idx]

    if '-t' in sys.argv:
        idx = sys.argv.index('-t')
        tamofile = sys.argv[idx+1]
        while len(sys.argv) > idx + 2 and  sys.argv[idx+2].isdigit():
            tamoidxs.append(int(sys.argv[idx+2]))
            del sys.argv[idx+2]
        del sys.argv[idx]
        del sys.argv[idx]

    if '-f' in sys.argv:
        idx = sys.argv.index('-f')
        factor = float(sys.argv[idx+1])
        del sys.argv[idx]
        del sys.argv[idx]

    if '-X' in sys.argv:
        idx = sys.argv.index('-X')
        del sys.argv[idx]
        find_best = 1

    if '-TF' in sys.argv:
        idx = sys.argv.index('-TF')
        del sys.argv[idx]
        TRANSFAC = 1

    if '-V' in sys.argv:
        idx = sys.argv.index('-V')
        del sys.argv[idx]
        verbose = 1

    if '-sub' in sys.argv:
        idx = sys.argv.index('-sub')
        subfsa = 'YES'
        del sys.argv[idx]

    if '-MAT' in sys.argv:  #Not Used
        idx = sys.argv.index('-MAT')
        print "## -MAT not available"
        TFMAT = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    if '-bysite' in sys.argv:
        idx = sys.argv.index('-bysite')
        method = 'bysite'
        del sys.argv[idx]
        if sys.argv[idx].isdigit():
            bysitek = int(sys.argv[idx])
            del sys.argv[idx]

    if '-bysitef' in sys.argv:
        idx = sys.argv.index('-bysitef')
        method = 'byfact'
        del sys.argv[idx]
        if sys.argv[idx].isdigit():
            byfactk = int(sys.argv[idx])
            del sys.argv[idx]

    for _method in ['spec', 'chi2', 'E_seq', 'MNCP', 'ROC_AUC',
                    'GO', 'feature', 'dist', 'binomial', 'overrep']:
        _flag = '-%s'%_method
        if _flag in sys.argv:
            idx = sys.argv.index(_flag)
            method = _method
            del sys.argv[idx]

    theProbeSet = ProbeSet(genome)
    theProbeSet.factor = factor
    if reducefile:
        ids = theProbeSet.ids_from_file(reducefile)
        theProbeSet.reduce_probelist(ids)


    NMfile = sys.argv[1]

    if subfsa:
        subD = Fasta.load(NMfile)
        print '# Substituting'
        for id,seq in subD.items():
            theProbeSet.probes[id] = seq

    if   N_file:
        if N_file < 11:
            for N in range(11,N_file-1,-1):
                motifs.extend(top_nmers_fsa(N,NMfile))
        else:
            for N in range(N_file,7,-1):
                t = top_nmers_fsa(N,NMfile)
                motifs.extend(top_nmers_fsa(N,NMfile))
    elif N_genome:
        if N_genome < 11:
            for N in range(11,N_genome-1,-1):
                motifs.extend(top_nmers_seqs(N,theProbeSet.seqs_from_idfile(NMfile,'quiet')))
        else:
            for N in range(N_genome,7,-1):
                motifs.extend(top_nmers_seqs(N,theProbeSet.seqs_from_idfile(NMfile,'quiet')))

    else:
        motifs.extend(sys.argv[2:])

    #Por que?
    #for i in range(len(motifs)-1,-1,-1):
    #    if not motifs[i].isalpha():
    #        del motifs[i]

    if toMotif:
        _t = motifs[:]
        motifs = []
        for motif in _t:
            m = MotifTools.Motif()
            m.compute_from_text(motif,0.1)
            motifs.append(m)
    elif tamofile:
        if 0 and re.search('\.ace$',tamofile):
            from TAMO import AlignAce
            motifs = AlignAce.AlignAce(tamofile).motifs
        else:
            tamomotifs = MotifTools.load(tamofile)
        for tamoidx in tamoidxs:
            motifs.append(tamomotifs[tamoidx])
        if not tamoidxs:
            motifs.extend(tamomotifs)
    if TFMAT: #Not used
        if not ('metaTF' in dir()):
            import metaTF
        #motifs.append(metaTF.matid2newmotif(TFMAT,bg={'A':.228, 'C':.272, 'G':.272, 'T':.228 }))
        motifs.append(metaTF.matid2newmotif(TFMAT,bg={'A':.31, 'C':.19, 'G':.19, 'T':.31 }))

    if TRANSFAC: #Not used
        tf_motifD = tfmats()
        for key,m in tf_motifD.items():
            if 1 or key.find('human') >= 0: #Remove '1 or' to limit to human
                m.source = key
                motifs.append(m)

    print "#Searching %d motifs"%len(motifs)
    ids = Fasta.keys(NMfile)
    if find_best:
        for motif in motifs:
            (p,frac) = theProbeSet.best_p_value(motif,ids,'verbose')
            sys.stdout.flush()
            seqsT = motif.bestseqs(frac * motif.maxscore)
            bestseqs = map(lambda d:d[1], seqsT)
            newMotif = MotifTools.Motif(bestseqs,motif.background)
            print 'Best: %-30s %e %5.4f'%(motif,p,frac)
            print newMotif
            newMotif._print_p()
            print 'New:'
            newMotif._print_bits()
            print 'Orig:'
            motif._print_bits()
    else:
        scores = []
        i = -1
        for motif in motifs:
            if   method == 'E':      p = theProbeSet.p_value(motif,ids,         'verbose',omitlist=omits)
            elif method == 'spec':   p = theProbeSet.church (motif,ids,         'verbose')
            elif method == 'bysite': p = theProbeSet.E_site (motif,ids,bysitek, 'verbose')
            elif method == 'byfact': p = theProbeSet.E_sitef(motif,ids,byfactk, 'verbose')
            elif method == 'overrep':p = theProbeSet.binomial(motif,ids,        'verbose')
            elif method == 'binomial':p= theProbeSet.binomial(motif,ids,        'verbose')
            elif method == 'chi2':   p = theProbeSet.E_chi2 (motif,ids,None,    'verbose')
            elif method == 'byseq':  p = theProbeSet.E_seq  (motif,ids,         'verbose')
            elif method == 'MNCP':   p = theProbeSet.MNCP   (motif,ids,         'verbose')
            elif method == 'ROC_AUC':p = theProbeSet.ROC_AUC(motif,ids,         'verbose')
            elif method == 'GO':     p = theProbeSet.GO     (motif,ids,         'verbose')
            elif method == 'feature':p = theProbeSet.feature_check(motif,ids,'vebose')
            elif method == 'dist':   p = theProbeSet.distribution(motif,ids)
            else:
                print "Unknown Method: %s"%(method)
                sys.exit(1)
            if verbose and type(motif) != type('') and p<0.1:
                #p = theProbeSet.p_value(motif,ids,'verbose') #That's right: again!
                i = i + 1
                motif.church = theProbeSet.church(motif,ids,'verbose',-1)
                motif.pvalue = p
                motif.nseqs  = theProbeSet._count_matching_probes_motif_fast(motif, ids, factor)
                MotifTools.print_motif(motif,20,i)
            scores.append((p,motif))
        if gapped:
            scores.sort()
            print "Top 10 Motifs"
            for i in range(min(10,len(scores))):
                score, motif = scores[i]
                if score == 1: continue
                print "Motif %2d  w: %2d   p: %4.2e  %s"%(i,len(motif),score,motif)
    sys.exit(0) #Avoid ridiculous cleanup times
            

class ProbeSet:
    '''
    class Probeset -- Class that represents groups of discrete genome sequences (such as may be represented
                      on a microarray, or the upstream regions of all genes) in such a way is to facilitate
                      computing motif scores and metrics.  See module documentation and documentation of
                      member functions.
    '''
    def __init__(self,genome='YEAST',factor=0.7):
        """
        ProbeSet.__init__(genome='YEAST',factor=0.7) -- 'genome' can be "YEAST" or a fasta-formatted file.
                                                        'factor' specifies the PSSM threshold in terms of
                                                        its maximum possible value.
        """
        self.factor = factor
        if genome=='YEAST':
            self.fastafile    = TAMO.paths.Whiteheaddir + 'Yeast6kArray/yeast_Young_6k.fsa'
            TAMO.paths.CHECK(self.fastafile,'Whitehead')
            self.xlatefile    = ''
            self._loadfiles()
        else:
            self.fastafile    = genome
            self.xlatefile    = ''
            self._loadfiles()
        if 0:
            print 'No Defaults for genome ',genome
        self.c_Probelists_str = {}
        self.totalbases = len(''.join(self.probes.values()))

    def _loadfiles(self):
        """
        ProbeSet._loadfiles() -- [Utility function]  Load the sequence file or pickled ProbeSet object.
        """
        self.probes       = {}
        self.gene2NM      = {}
        if   os.path.exists(self.fastafile + '.pickle'):
            FID = open(self.fastafile + '.pickle','r')
            self.probes = pickle.load(FID)
        elif os.path.exists(self.fastafile + '.shelf'):
            self.probes = shelve.open(self.fastafile + '.shelf','r')
            _t = {}
            for k in self.probes.keys():
                _t[k] = self.probes[k]
            self.probes = _t
        else:
            self.probes = Fasta.load(self.fastafile)

        if not  os.path.exists(self.fastafile + '.pickle'):
            print 'Writing %s...'%(self.fastafile + '.pickle')
            FID = open(self.fastafile + '.pickle','w')
            pickle.dump(self.probes,FID)
        
        if self.xlatefile:
            if os.path.exists(self.xlatefile + '.shelf'):
                self.gene2NM = shelve.open(self.xlatefile + '.shelf','r')
            else:
                self.gene2NM = shelve.open(self.xlatefile + '.shelf','r')
                FID = open(self.xlatefile,'r')
                for line in FID.readlines():
                    gene,NM = line.split()[0:2]
                    self.gene2NM[gene] = NM
                FID.close()

        for key in self.probes.keys():
            self.probes[key]=self.probes[key].upper()
            
    def _count_matching_probes_motif_fast(self, motif, NMlist=[], factor=None):
        """
        ProbeSet._count_matching_probes_motif_fast(motif, NMlist=[], factor=None) -- [Utility, C++ Interface]
                 Count matching probes.  Optional NMlist limits search to entries with names in list, factor
                 specifies PSSM threshold.
        """
        c_PSSM = MDsupport.SeqMat(motif.width)
        LjT = zip(['A','C','G','T'],[0,1,2,3])
        for i in range(motif.width):
            for L,j in LjT:
                c_PSSM.set(i,j, motif.ll[i][L])
        if NMlist == []: NMlist = self.probes.keys()
        count = 0
        #print "Scanning with",motif
        if factor != None:
            threshold = motif.maxscore * factor
        elif motif.threshold == None:
            threshold = motif.maxscore * self.factor
        else:
            threshold = motif.threshold
        i = 0
        for key in NMlist:
            i = i + 1
            seq   = self.probes[key]
            try:
                bestscore = c_PSSM.scanbest(seq)
            except:
                print key, self.probes[key]
                sys.exit(1)
            if bestscore >= threshold:
                count = count + 1
        return(count)

    def _count_matching_sites_motif_fast(self, motif, ceil=1, NMlist=[], factor=None):
        """
        ProbeSet._count_matching_sites_motif_fast( motif, ceil=1, NMlist=[], factor=None) -- [Utility, C++ interface]
                 Tally # of maching sites.  ceil specifies maximum number of sites to tally per probe/seqeunce entry
                 Optional NMlist limits search to entries with names in list, factor speciies PSSM threshold in
                 terms of motif maximum possible score.
        """
        c_PSSM = MDsupport.Motif2c_PSSM(motif)
        if NMlist == []: NMlist = self.probes.keys()
        if factor != None:
            threshold = motif.maxscore * factor
        elif motif.threshold == None:
            threshold = motif.maxscore * self.factor
        else:  threshold = motif.threshold
        starts = []
        count = 0.0
        for key in NMlist:
            seq   = self.probes[key]
            hits  = c_PSSM.matchstarts(seq,threshold)
            if len(hits) >= ceil: count = count + 1
            #count = count + min(ceil,len(hits))/float(ceil)
        return(int(count))

    def _countD_matching_sites_motif_fast(self, motif, NMlist=[], factor=None, Nooverlap=None):
        """
        ProbeSet._countD_matching_sites_motif_fast( motif, NMlist=[], factor=None, Nooverlap=None) -- [Utility, C++]
                 Return a dictionary of motif counts per probe, such as might be used for a histogram or for
                 computing Chi-square.
                 Optional NMlist limits search to entries with names in list, factor speciies PSSM threshold in
                 terms of motif maximum possible score.  Nooverlap flag prohibits counting of sites within a
                 motif's width of one another.
        """
        #Nooverlap = None
        D = {}
        c_PSSM = MDsupport.Motif2c_PSSM(motif)
        if NMlist == []: NMlist = self.probes.keys()
        if factor != None:
            threshold = motif.maxscore * factor
        elif motif.threshold == None:
            threshold = motif.maxscore * self.factor
        else:  threshold = motif.threshold
        starts = []
        count = 0.0
        w = motif.width
        for key in NMlist:
            seq   = self.probes[key]
            hits  = c_PSSM.matchstarts(seq,threshold)
            if Nooverlap and len(hits) > 1:
                ans = [hits[0]]
                for i in range(1,len(hits)):
                    if hits[i] >= ans[-1] + w: ans.append(hits[i])
                    #else: print '(%d - %d) = %d < %d'%(hits[i],w,(hits[i]-w),hits[i-1])
                hits = ans
            count = len(hits)
            #if count: print key,count
            try:             D[count] = D[count] + 1
            except KeyError: D[count] = 1
        return(D)

    def _probe_scores(self, motif, NMlist=[],base_threshold=-100):
        """
        ProbeSet._probe_scores( motif, NMlist=[],base_threshold=-100) -- [Utility, C++]
                 Return a list of (score, seq_id) tuples, where score is the score of the best
                 matching site on the sequence with id "seq_id" in the ProbeSet.
        """
        c_PSSM = MDsupport.SeqMat(motif.width)
        LjT = zip(['A','C','G','T'],[0,1,2,3])
        for i in range(motif.width):
            for L,j in LjT:
                c_PSSM.set(i,j, motif.ll[i][L])
        if NMlist == []: NMlist = self.probes.keys()
        ansT = []
        count = 0
        #print "Scanning with",motif
        #threshold = motif.maxscore * self.factor
        threshold = base_threshold
        i = 0
        for key in NMlist:
            i = i + 1
            seq   = self.probes[key]
            bestscore = c_PSSM.scanbest(seq)
            if bestscore >= threshold:
                count = count + 1
                ansT.append((bestscore,key))
                #print count,key,bestscore
            if (i%500) == 0:
                #print '.',;sys.stdout.flush()
                pass
        return(ansT)

    def _count_matching_sites(self, motif, threshold):
        """
        ProbeSet._count_matching_sites( motif, threshold) -- [Utility, C++, Deprecated]
                 Use _count_matching_sites_motif_fast instead.  See it's docstring.
                 
        """
        c_PSSM = MDsupport.Motif2c_PSSM(motif)
        starts = []
        for seq in seqs:
            starts.extend(c_PSSM.matchstarts(seq,threshold))
        return len(starts)

    def _count_matching_probes_motif(self, motif, NMlist=[]):
        """
        ProbeSet._count_matching_probes_motif( motif, NMlist=[]) -- [Utility, C++, Deprecated]
                 Use _count_matching_probes_motif_fast instead.  See it's docstring.
        """
        if NMlist == []: NMlist = self.probes.keys()
        count = 0
        print "# Scanning with",motif
        threshold = motif.maxscore * self.factor
        #print threshold, motif.maxscore
        #motif._print_ll()
        i = 0
        for key in NMlist:
            i = i + 1
            seq   = self.probes[key]
            (matches,endpoints,scores) = motif._scan(seq,threshold)
            if (i%500) == 0:
                #print '.',;sys.stdout.flush()
                pass
            if len(matches) > 0:
                count = count + 1
                #print count,key,zip(matches,scores)
        print
        return(count)

    def _count_matching_probes_regex_fast(self, ambig, NMlist=[]):
        """
        ProbeSet._count_matching_probes_regex_fast( ambig, NMlist=[]) -- [Utility, C++]
                 Count how many probes contain a substring that matches the regular expression. 
        """
        #key = '.'.join(NMlist)
        key = len(NMlist)
        try: key = key + NMlist[0] + NMlist[-1]
        except: pass
        if self.c_Probelists_str.has_key(key):
            P = self.c_Probelists_str[key]
        else:
            P = MDsupport.Probelist_str()
            for NM in NMlist:
                P.append(self.probes[NM])
            self.c_Probelists_str[key] = P
        count = P.count_re_matches(ambig)
        return count

    def _count_matching_probes_regex(self, ambig, NMlist=[]):
        """
        ProbeSet._count_matching_probes_regex( ambig, NMlist=[]) -- [Utility, C++, deprecated]
                 Use _count_matching_probes_regex_fast instead.  See its docstring.
        """
        if NMlist == []: NMlist = self.probes.keys()
        count = 0
        site_re = re.compile(AmbigToRegExp(ambig))
        for key in NMlist:
            seq   = self.probes[key]
            seqrc = revcomp(seq)
            if (site_re.search(seq) or site_re.search(seqrc)):
                #print 'Match: ',key,site_re.search(seq),site_re.search(seqrc)
                count = count + 1
        return(count)

    def count_matching_sites(self, ambig, ceil=1,NMlist=[], thresh=None):
        """
        ProbeSet.count_matching_sites( ambig, ceil=1,NMlist=[], thresh=None) -- [Utility, Dispatch]
                 Confirm that ambig is a motif object (as opposed to a regular expression) and invoke
                 _count_matching_sites_motif_fast().   See its docstring.
        
        """
        if NMlist == []: NMlist = self.probes.keys()
        if type(ambig) == type("STRING"):
            raise "MotifReqd4MatchingSites"
        elif (ambig.__class__ == MotifTools.Motif) or type(ambig) == type(MotifTools.Motif()):
            count = self._count_matching_sites_motif_fast(ambig, NMlist, ceil, thresh)
        else:
            print ambig.__class__
            print ambig
            raise "UnknownMotiforAmbiguitySequence"

        return(count)

    def countD_matching_sites(self, ambig, NMlist=[], thresh=None,Nooverlap=None):
        """
        ProbeSet.countD_matching_sites( ambig, NMlist=[], thresh=None,Nooverlap=None) -- [Utility, Dispatch]
                 Confirm that ambig is a motif object (as opposed to a regular expression) and invoke
                 _countD_matching_sites_motif_fast.  See its documentation.
        """
        if NMlist == []: NMlist = self.probes.keys()
        if type(ambig) == type("STRING"):
            raise "MotifReqd4MatchingSites"
        elif (ambig.__class__ == MotifTools.Motif) or type(ambig) == type(MotifTools.Motif()):
            count = self._countD_matching_sites_motif_fast(ambig, NMlist, thresh,Nooverlap)
        else:
            print ambig.__class__
            print ambig
            raise "UnknownMotiforAmbiguitySequence"

        return(count)

    def count_matching_probes(self, ambig, NMlist=[], thresh=None): 
        """
        ProbeSet.count_matching_probes( ambig, NMlist=[], thresh=None): --  [Utility, Dispatch]
                 Check whether ambig is a Motif object or a regular expression, and invoke the
                 appropriate routine to use ambig to count the number of matching probes. 
        """
        if NMlist == []: NMlist = self.probes.keys()
        if type(ambig) == type("STRING"):
            count = self._count_matching_probes_regex(ambig, NMlist)
        elif (ambig.__class__ == MotifTools.Motif) or type(ambig) == type(MotifTools.Motif()):
            count = self._count_matching_probes_motif_fast(ambig, NMlist, thresh)
            #count = self._count_matching_sites_motif_fast(ambig, NMlist, thresh)
        else:
            print ambig.__class__
            print ambig
            raise "UnknownMotiforAmbiguitySequence"

        return(count)

    def Enrichment(self,ambig, NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.Enrichment(ambig, NMlist,verbose='',factor=None,omitlist=[]) --
                 Computes Enrichment score as described in Harbison & Gordon et al. Nature 2004
                 Optional factor is expressed in terms of fraction of the maxiumum possible
                 PSSM ("ambig") score.  Result is computed by using the hypergeometric distribution
                 to compute the probability of observing the number entries in NMlist that
                 match the motif by considering them has drawn randomly from the ProbeSet without
                 replacement.
        """
        return self.p_value(ambig, NMlist,verbose   ,factor     ,omitlist)
    def E      (self,ambig, NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.E      (ambig, NMlist,verbose='',factor=None,omitlist=[]) --
                 Synonym for ProbeSet.Enrichment.  See its docstring
        """
        return self.p_value(ambig, NMlist,verbose   ,factor     ,omitlist)
    def p_value(self,ambig, NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.p_value(ambig, NMlist,verbose='',factor=None,omitlist=[]) -- 
                 Synonym for ProbeSet.Enrichment.  See its docstring
        """
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        
        total          = len(all)
        numpicked      = len(subset)
        if not subset: numfound = 0
        else:          numfound = self.count_matching_probes(ambig,subset,factor)
        if numfound:
            numinteresting = self.count_matching_probes(ambig,all,factor)
            p = Arith.hypgeomsummore(numinteresting,total,numpicked,numfound)
        else:
            numinteresting = -1
            p = 1
        if verbose:
            print "E:     %4d / %4d, %3d / %4d found.       p: %4.2e  %s"%(
                numinteresting,total, numfound,numpicked,p,ambig),
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print
        return(p)

    def p_value_noomit(self,ambig, NMlist,verbose='',factor=None):
        """
        ProbeSet.p_value_noomit(ambig, NMlist,verbose='',factor=None) -- [Deprecated]
                 Synonym for ProbeSet.Enrichment.  See its docstring
        """
        total          = len(self.probes.keys())
        numpicked      = len(NMlist)
        numfound       = self.count_matching_probes(ambig,NMlist,factor)
        if numfound:
            numinteresting = self.count_matching_probes(ambig,[],factor)
            p = Arith.hypgeomsummore(numinteresting,total,numpicked,numfound)
        else:
            numinteresting = -1
            p = 1
        if verbose:
            print "E:  %4d / %4d, %3d / %4d found.  p: %4.2e  %s"%(
                numinteresting,total, numfound,numpicked,p,ambig),
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print
        return(p)


    def best_Enrichment(self,ambig, NMlist,verbose=''):
        """
        ProbeSet.best_Enrichment(ambig, NMlist,verbose='') -- 
                 Find the scoring threshold at which the motif ("ambig") has the optimal Enrichment score
                 (see ProbeSet.Enrichment).  Returns (best_enrichment_score, best_factor) tuple.
        """
        return best_p_value(self,ambig, NMlist,verbose)
    def best_pvalue(self,ambig, NMlist,verbose=''):
        """
        ProbeSet.best_pvalue(ambig, NMlist,verbose='') --
                 Find the scoring threshold at which the motif ("ambig") has the optimal Enrichment score
                 (see ProbeSet.Enrichment)  Returns (best_enrichment_score, best_factor) tuple.  (Synonym
                 for best_Enrichment)
        """
        return best_p_value(self,ambig, NMlist,verbose)

    def best_p_value(self,ambig, NMlist,verbose=''):
        """
        ProbeSet.best_p_value(ambig, NMlist,verbose='') --
                 Synonym for best_Enrichment.  See its docstring. 
        """
        maxprecision = 2
        best = 0.5
        bestp = 1
        _factor = self.factor
        for precision in range(1,maxprecision+1):
            f = math.pow(10,-precision)
            for factor in map(lambda x,b=best,f=f:b + f*x, range(-5,5)):
                self.factor = factor
                total          = len(self.probes.keys())
                numpicked      = len(NMlist)
                numfound       = self._count_matching_probes_motif_fast(ambig, NMlist, factor)
                if numfound:
                    numinteresting = self._count_matching_probes_motif_fast(ambig)
                    p = Arith.hypgeomsummore(numinteresting,total,numpicked,numfound)
                else:
                    numinteresting = -1
                    p = 1
                if verbose:
                    print "%5.3f  Of %4d / %4d, %3d / %4d found.  p: %4.2e  %s"%(
                        factor, numinteresting,total, numfound,numpicked,p,ambig)
                    sys.stdout.flush()
                if p < bestp:
                    best  = factor
                    bestp = p
        self.factor = _factor
        return(bestp,best)

    def church(self,ambig,preNMlist,verbose='',topcount=100,omitlist=[]):
        """
        ProbeSet.church(ambig,preNMlist,verbose='',topcount=None,omitlist=[]) --
                 The "Group Specificity Score" as used by AlignACE (Hughes et al. JMB 2000)
                 Optional topcount is used as "100" in AlignACE.  A negative fractional value for
                 topcount is used to specify the count in terms of the number of sequences
                 in NMlist.
        """
        if omitlist:
            NMlist = [x for x in preNMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            NMlist = preNMlist
            all    = self.probes.keys()
        scoresT = self._probe_scores(ambig,all)
        scoresT.sort()
        scoresT.reverse()

        if topcount<0:      topcount = int(-topcount*len(NMlist))

        #Figure out whether to take top 100 or more, as in a case of a tie between 100 & 101, etc...
        if len(scoresT) > topcount:
            for i in range(topcount,len(scoresT)):
                if not scoresT[i][0] == scoresT[topcount][0]: break
            scoresT = scoresT[0:i]

        #Cutoff at threshold
        if 0:  #No threshold in AlignACE?
            threshold = ambig.maxscore * self.factor
            last = 0
            for score,id in scoresT:
                if score > threshold:
                    last = last + 1
                else:
                    break
            scoresT = scoresT[0:last]

        total     = len(all)   #Total number of probes
        numinteresting = len(NMlist)
        numpicked  = len(scoresT)              #Number of genome matches (particulary if > 100)

        numfound = 0
        for i in range(len(scoresT)):
            if scoresT[i][1] in NMlist: numfound = numfound + 1

        if numfound:
            p = Arith.hypgeomsummore(numinteresting,total,numpicked,numfound)
        else:
            p = 1
        if verbose:
            print "Church %4d / %4d, %3d / %4d top-ranked.  p: %4.2e  %s"%(
                numinteresting,total, numfound,numpicked,p,ambig),
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print
        return(p)
        
    def E_site(self,ambig, NMlist,ceil=1,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.E_site(ambig, NMlist,ceil=1,verbose='',factor=None,omitlist=[]) --
                 Like the enrichment score, but computes the hypergrometric on the number of sites.
                 Not recommended, because the "total" and "numpicked" are based on number of probes,
                 so the metric intrinsically contains an apple vs. orange comparison.
        """
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        
        total          = len(all)
        numpicked      = len(subset)
        numfound       = self.count_matching_sites(ambig,subset,ceil,factor)
        #Pnumfound      = self.count_matching_probes(ambig,subset,factor)
        if numfound:
            numinteresting = self.count_matching_sites(ambig,all,ceil,factor)
            #Pnuminteresting = self.count_matching_probes(ambig,all,factor)
            p = Arith.hypgeomsummore(numinteresting,total,numpicked,numfound)
        else:
            numinteresting = -1
            p = 1
        if verbose:
            print "Es: Of %4d / %4d, %3d / %4d found.       p: %4.2e  k=%d %s"%(
                numinteresting,total, numfound,numpicked,p,ceil,ambig),
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print
        return(p)

    def E_chi2(self,ambig, NMlist,bins=None,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.E_chi2(ambig, NMlist,bins=None,verbose='',factor=None,omitlist=[]) --
                 Metric that compares the distribution of number of probes with multiple occurrences
                 as compared to the distribution seen the the genome at large.  Not completely stable,
                 because counts are often too small for reliable chi2.
        """
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        if not bins: bins = [0,1,2,3,4,5]

        #Facilitate placing counts into bins with index lookup
        val2binD = {}
        for i in range(10):
            if i >= max(bins):
                val2binD[i]=len(bins)-1 #last bin
            else:
                for j in range(len(bins)-1):
                    if (i >= bins[j]) and (i<bins[j+1]):
                        val2binD[i]=j
                        continue
        #for val,bin in val2binD.items():
        #    print val, bin

        #Measure actual occurences, first in bound probes, then in genome
        total          = len(all)
        numpicked      = len(subset)
        boundD         = self.countD_matching_sites(ambig,subset,factor)
        if boundD:
            genomeD    = self.countD_matching_sites(ambig,all,factor)
            #Initialize arrays to hold binned values
            #print boundD
            #print genomeD
            boundbins  = [0] * len(bins)
            genomebins = [0] * len(bins)
            boundtot   = 0.0
            genometot  = 0.0
            #Place counts into bins
            for BIN,D in [(boundbins,boundD), (genomebins,genomeD)]:
                for value,count in D.items():
                    try:    idx = val2binD[value]
                    except: idx = len(bins)-1
                    BIN[idx] = BIN[idx] + count
            for count in boundbins:  boundtot  = boundtot  + count
            for count in genomebins: genometot = genometot + count
            expectbins = [float(x)/genometot*boundtot for x in genomebins]
            for i in range(len(bins)-1,1,-1):
                if (expectbins[i] < 5) or (boundbins[i] < 5):
                    genomebins[i-1] = genomebins[i-1] + genomebins[i]
                    expectbins[i-1] = expectbins[i-1] + expectbins[i]
                    boundbins [i-1] = boundbins [i-1] +  boundbins[i]
                    del expectbins[i], boundbins[i], genomebins[i]
            #Compute Chi2
            if verbose:
                print "[ %s ]:: [ %s ] vs. [ %s ] = "%(
                    ', '.join(['%d'%x for x in bins]),
                    ', '.join(['%d'%x for x in genomebins]),
                    ', '.join(['%d'%x for x in boundbins])),

            chisq,pvalue=stats.lchisquare(boundbins,expectbins)
            #chisq,pvalue=stats.lchisquare(genomebins,boundbins)
            p = pvalue
        else:
            p = 1
        if verbose:
            try:
                print "%4.3e %8.2f"%(p,chisq),
            except:
                print p,chisq,len(genomebins)-1
                print stats.lchisqprob(chisq,len(genomebins)-1)
                sys.exit(1)
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print ambig.oneletter
        return(p)

                
    def distribution(self,motif,NMlist,):
        """
        ProbeSet.distribution(motif,NMlist,) --
                 Prints the "id   best_score, bestscore/maxscore" information for all probes in NMlist,
                 followed by some general statistical information.  If Stats.stats is installed,
                 also attempts KS test of distribution of scores within bound probes versus the genome. 
        """
        c_PSSM = MDsupport.Motif2c_PSSM(motif)
        bests, ratios, allratios = [], [], []
        for id in NMlist:
            #bestscore = c_PSSM.scanbest(self.probes[id])
            bestscore,seq = motif.bestscanseq(self.probes[id])
            rcseq = MotifTools.revcomplement(seq)
            if motif.score(seq,'FWD') < motif.score(rcseq,'FWD'):
                seq = rcseq
            ratio     = bestscore/motif.maxscore
            print '@ %-15s  %6.2f %6.3f %s %s'%(id,bestscore,ratio,seq,motif.oneletter)
            bests.append(bestscore)
            ratios.append(ratio)
        for id in self.probes.keys():
            bestscore = c_PSSM.scanbest(self.probes[id])
            ratio     = bestscore/motif.maxscore
            allratios.append(ratio)
        P = 1
        try:
            W,P = stats.ks_2samp(allratios,ratios)
        except:
            W   = 1
        bave,bstd = Arith.avestd(bests)
        rave,rstd = Arith.avestd(ratios)
        print '-'*40
        print '%-15s %6.2f +/- %3.1f    %6.3f +/- %5.3f  %5.3f %5.3e %s'%(
            'ALL', bave, bstd, rave, rstd, W, P, motif.oneletter )
        print '-'*40,'\n'
        return P  #Parent sometimes expects a pvalue

    def matching_ids(self,ambig,NMlist=[],factor=None):
        """
        ProbeSet.matching_ids(ambig,NMlist,factor=None) --
                 Return a list of sequence ids among NMlist (probe or gene ids) that match the motif.
                 Optional factor is expressed in terms of fraction of the maxiumum possible PSSM score.
        """
        matches = []
        if not NMlist: NMlist = self.probes.keys()
        if type(ambig) == type("STRING"):
            pass
            for id in NMlist:
                if re.search(ambig,self.probes[id]):
                    matches.append(id)
        elif (ambig.__class__ == MotifTools.Motif) or type(ambig) == type(MotifTools.Motif()):
            motif = ambig
            if factor != None:            threshold = motif.maxscore * factor
            elif motif.threshold == None: threshold = motif.maxscore * self.factor
            else:                         threshold = motif.threshold
            c_PSSM = MDsupport.Motif2c_PSSM(motif)
            for id in NMlist:
                bestscore = c_PSSM.scanbest(self.probes[id])
                if bestscore >= threshold:
                    matches.append(id)
        return matches

    def matches_with_features(self,ambig,NMlist,featurelist,featurename='fE',verbose='',factor=None):
        """
        ProbeSet.matches_with_features(ambig,NMlist,featurelist,featurename='fE',verbose='',factor=None) --
                 Check Enrichment of motif among ids in featurelist that have something "interesting" about
                 them.  For example, you would use this function if you wanted to ask if a motif has a good
                 Enrichment score among probes identified in a microarray experiment that are, say,
                 telomeric, or perhaps known to be important during cell-cycle.
                 NMlist and featurelist are both lists of gene/probe ids.  Featurename is used only for
                 visual reference when printing the results.  Optional factor is expressed in terms
                 of fraction of the maxiumum possible PSSM score.
        """
        all        = featurelist
        subset     = [x for x in NMlist if x in featurelist]
        matches    = self.matching_ids(ambig,featurelist,factor)
        submatches = [x for x in matches if x in subset]
        if submatches:
            p = Arith.hypgeomsummore(len(matches), len(all),
                                     len(subset),  len(submatches))
            pG = Arith.hypgeomsummore(len(all),     len(self.probes.keys()),
                                     len(NMlist),  len(subset))
                                                       
        else:
            p = 1
            pG = 1
        if verbose:
            if p/pG < 0.1:
                print "%-29s  %5d / %5d, %5d / %5d found.  ( %4.2e / %4.2e ) = p/pG: %4.2e %s"%(
                    featurename+':',
                    len(matches),len(all),len(submatches),len(subset),p,pG,p/pG,ambig),
                if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                    print ambig.source
                else:
                    print
        return(p)
                    
    def GO(self,ambig,NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.GO(ambig,NMlist,verbose='',factor=None,omitlist=[]) --
                 Check whether the subset of NMlist ids that match the motif have significant
                 categorization according to the Gene Ontology (GO) database.
        """
        if 'GO' not in dir():
            from TAMO.DataSources import GO
        if omitlist: subset = [x for x in NMlist if not x in omitlist]
        else:        subset = NMlist
        if omitlist: all    = [x for x in self.probes.keys() if not x in omitlist]
        else:        all    = self.probes.keys()
        
        matches = self.matching_ids(ambig,all,factor)
        results = GO.probelist2categories(matches,1)
        if not results: result = ''
        else          : result = results[0]
        print "GO: Of %3d matching probes, %s %s"%(len(matches),result,ambig)

    def MNCP(self,ambig, NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.MNCP(ambig, NMlist,verbose='',factor=None,omitlist=[]) --
                 Compute MNCP as decribed in Clarke ND, Granek JA  Bioinformatics 2003
        """
        if 'calc_MNCP' not in dir():
            from TAMO.util.Clarke.narke_test import narke_score as calc_MNCP
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        motif = ambig
        if factor != None:            threshold = motif.maxscore * factor
        elif motif.threshold == None: threshold = motif.maxscore * self.factor
        else:                         threshold = motif.threshold
        c_PSSM = MDsupport.Motif2c_PSSM(motif)
        D = {}
        w = motif.width
        for id in all:
            #D[id] = c_PSSM.sumscoresabove(self.probes[id],threshold)
            starts = c_PSSM.matchstarts(self.probes[id],threshold)
            if starts:
                ans = [starts[0]]
                for start in starts[1:]:
                    if start >= ans[-1] + w: ans.append(start)
                D[id] = len(ans)
            else: D[id] = 0
        subsetscores = [D[x] for x in subset]
        allscores    = D.values()
        nmcp = calc_MNCP(subsetscores,allscores,1)
        if verbose:
            print "MNCP:  %-20s mn: %8.4f "%(motif,nmcp)
        return nmcp

    def ROC_AUC(self,ambig, NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.ROC_AUC(ambig, NMlist,verbose='',factor=None,omitlist=[]) -- 
                 Compute ROC_AUC as decribed in Clarke ND, Granek JA  Bioinformatics 2003
        """
        if 'ROC_auc' not in dir():
            from TAMO.util.Clarke.ROC_AUC    import ROC_auc, normalize_auc
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        motif = ambig
        if factor != None:            threshold = motif.maxscore * factor
        elif motif.threshold == None: threshold = motif.maxscore * self.factor
        else:                         threshold = motif.threshold
        c_PSSM = MDsupport.Motif2c_PSSM(motif)
        w = motif.width
        D = {}
        for id in all:
            #D[id] = c_PSSM.sumscoresabove(self.probes[id],threshold)
            starts = c_PSSM.matchstarts(self.probes[id],threshold)
            if starts:
                ans = [starts[0]]
                for start in starts[1:]:
                    if start >= ans[-1] + w: ans.append(start)
                D[id] = len(ans)
            else: D[id] = 0
        subsetscores = [D[x] for x in subset]
        otherscores    = [D[x] for x in all if x not in subset]
        auc = ROC_auc(subsetscores,otherscores,1)
        nauc = normalize_auc(auc,subsetscores,otherscores)
        if verbose:
            print "ROC_auc: %-20s ra: %8.4f %8.4f"%(motif,auc,nauc)
        return nauc
        
    def overrep(self,ambig,NMlist, verbose='',factor=None,omitlist=[],MAXVALUE=4):
        """
        ProbeSet.overrep(ambig,NMlist, verbose='',factor=None,omitlist=[],MAXVALUE=4) -- 
                 Computes Over-representation score.
                 Optional factor is expressed in terms of fraction of the maxiumum possible
                 PSSM ("ambig") score.

                 The score is computed by using the binomial distribution to compute the
                 probability of observing the number sites that match the motif among the sequences
                 in NMlist by considering them has drawn randomly from the ProbeSet with
                 replacement.  The expected fractions are based on the number of bases in the NMlist
                 probe/gene sequences and the total number of bases in the genome.

                 MAXVALUE is the maximum number of binding sites to count within a probe/gene
                 sequence.  A good value for 500-1000bp sequences is 4, but MAXVALUE should
                 be made very high if analyzing longer sequences. 
                 
        """
        return self.binomial(ambig,NMlist,verbose,factor,omitlist,MAXVALUE)
    
    def binomial(self,ambig,NMlist, verbose='',factor=None,omitlist=[],MAXVALUE=4):
        """
        ProbeSet.binomial(ambig,NMlist, verbose='',factor=None,omitlist=[],MAXVALUE=4) --
                 Synonym for ProbeSet.overrep.  See its docstring.
        """
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        #TOTAL
        total     = self.totalbases 
        numpicked = len(''.join([self.probes[x] for x in subset]))
        frac_exp  = float(numpicked)/float(total)

        if type(ambig) == type(''):
            ambig = MotifTools.Motif_from_text(ambig)
        
        countDsub = self.countD_matching_sites(ambig,subset,factor,'NoOverlap')
        if countDsub:
            countDall = self.countD_matching_sites(ambig,all,factor,'NoOverlap')
        
            numinteresting, numfound  = 0,0
            if countDall.has_key(0):
                del countDall[0] #We don't want to count "zero" occurrences
            if countDsub.has_key(0):
                del countDsub[0]
            for value,count in countDall.items(): numinteresting = numinteresting + min(value,MAXVALUE)*count
            for value,count in countDsub.items(): numfound       = numfound       + min(value,MAXVALUE)*count
            p = Arith.binomialsumtail(frac_exp,numinteresting,numfound)
        else:
            numinteresting = -1
            p = 1.0
        if verbose:
            print "Ov:  Observed %6d / %8d, expecting  %12.7f  ( %4.1f )  p: %4.2e %s"%(
                numfound,numinteresting,frac_exp,frac_exp*numinteresting,p,ambig),
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print
        return(p)
        

    def E_seq(self,ambig, NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.E_seq(ambig, NMlist,verbose='',factor=None,omitlist=[]) -- [Unstable]
                 Compute the Enrichment score, with these two changes:
                 The number of "draws" is the number of bases (both in the genome and
                 in the NMlist probe/gene sequences) and the motif counts are the number
                 of matching, but non-ovelapping binding sites.  In theory, this should
                 behave similarly to the overrep score with MAXVALUE set very high.
        """
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        #TOTAL
        total     = self.totalbases 
        numpicked = len(''.join([self.probes[x] for x in subset]))
        countDsub = self.countD_matching_sites(ambig,subset,factor,'NoOverlap')
        if countDsub:
            countDall = self.countD_matching_sites(ambig,all,factor,'NoOverlap')
            numinteresting, numfound  = 0,0
            if countDall.has_key(0):
                del countDall[0] #We don't want to count how many probe don't have the motif
            if countDsub.has_key(0):
                del countDsub[0]
            for value,count in countDall.items(): numinteresting = numinteresting + value*count
            for value,count in countDsub.items(): numfound       = numfound       + value*count
            p = Arith.hypgeomsummore(numinteresting,total,numpicked,numfound)
        else:
            numinteresting = -1
            p = 1
        if verbose:
            print "Eq:  %6d / %8d, %6d / %6d found.       p: %4.2e %s"%(
                numinteresting,total, numfound,numpicked,p,ambig),
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print
        return(p)
    def E_sitef(self,ambig, NMlist,k, verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.E_sitef(ambig, NMlist,k, verbose='',factor=None,omitlist=[]) --
                 Another variant of the Enrichment score, which requires at least
                 k occurences of a motif to occur on a probe for the probe to be
                 counted.  Each probe is counted k times. 
        """
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
            all    = [x for x in self.probes.keys() if not x in omitlist]
        else:
            subset = NMlist
            all    = self.probes.keys()
        #TOTAL
        total     = len(all)*k
        numpicked = len(subset)*k
        countDsub = self.countD_matching_sites(ambig,subset,factor,'NoOverlap')
        if countDsub:
            countDall = self.countD_matching_sites(ambig,all,factor,'NoOverlap')
            numinteresting, numfound  = 0,0
            if countDall.has_key(0):
                del countDall[0] #We don't want to count how many probe don't have the motif
            if countDsub.has_key(0):
                del countDsub[0]
            for value,count in countDall.items():
                if value > k: value = k
                numinteresting = numinteresting + value*count
            for value,count in countDsub.items():
                if value > k: value = k
                numfound = numfound + value*count
            p = Arith.hypgeomsummore(numinteresting,total,numpicked,numfound)
        else:
            numinteresting = -1
            p = 1
        if verbose:
            print "Esf:   %6d / %8d, %6d / %6d found.       p: %4.2e %s"%(
                numinteresting,total, numfound,numpicked,p,ambig),
            if type(ambig) != type('STRING') and type(ambig.source) == type('STRING'):
                print ambig.source
            else:
                print
        return(p)

    def frac(self,ambig, NMlist,verbose='',factor=None,omitlist=[]):
        """
        ProbeSet.frac(ambig, NMlist,verbose='',factor=None,omitlist=[]) --
                 Report and return the fraction of the probes/gene sequences in NMlist
                 that contain matches to the motif ("ambig") under the default
                 threshold score or the specified theshold according to factor (which
                 indicates the fraction of the maximum possible PSSM score to be used
                 as a threshold.
        """
        if omitlist:
            subset = [x for x in NMlist if not x in omitlist]
        else:
            subset = NMlist
        
        numpicked      = len(subset)
        numfound       = self.count_matching_probes(ambig,subset,factor)
        fraction = float(numfound)/float(numpicked)
        if verbose:
            print 'Frac: %4d / %4d f: %6.4f %s'%(
                numfound, numpicked, fraction, ambig)
        return fraction

        
    def ids_from_file(self,filename,QUIET=''):
        """
        ProbeSet.ids_from_file(filename,QUIET='') -- [Deprecated]
                 Fasta and txt file loader.  Use Fasta.keys() or Fasta.ids() instead.
        """
        ids = []
        fasta_flag = 0
        if filename.find('.fsa')>=0 or filename.find('.fasta')>=0:
            fasta_flag = 1
        FID = open(filename,'r')
        for line in FID.readlines():
            if fasta_flag and line[0] != '>':
                continue
            gene = line.split()[0].strip().replace('>','')
            if self.xlatefile and self.gene2NM.has_key(gene):
                gene = self.gene2NM[gene]
            if self.probes.has_key(gene):
                ids.append(gene)
            else:
                if not QUIET: print "# No Entry in genome for %s"%gene
        return(ids)

    def fsa_string_from_ids(self, ids, comments=[]):
        """
        ProbeSet.fsa_string_from_ids( ids, comments=[]) --
                 Prepare a text string from ids for
                 printing in Fasta format.  Similar functionality in Fasta module.
        """
        entries = []
        if type(ids) == type(''):
            ids = [ids]
        for i in range(len(ids)):
            s = ''
            id = ids[i]
            if not self.probes.has_key(id): continue
            comment = ''
            if i < len(comments): comment = comments[i]
            s = s + '>%s %s\n'%(id,comment)
            seq = self.probes[id]
            for i in range(0,len(seq),70):
                s = s + seq[i:i+70] + '\n'
            entries.append(s)
        s = ''.join(entries)
        return(s[:-1])

    def print_fsa_from_ids(self,ids):
        """
        ProbeSet.print_fsa_from_ids(ids) --
                 Print Fasta-formatted sequences to standard output.
                 Similar functionality also exists in the Fasta module.
        """
        for id in ids:
            print ">%s"%id
            seq = self.probes[id]
            for i in range(0,len(seq),70):
                print seq[i:i+70]


    def filter(self,probelist):
        """
        ProbeSet.filter(probelist) --
                 Return list of ids that are also represented in the ProbeSet object.
        """
        ans = []
        for probe in probelist:
            if self.probes.has_key(probe):
                ans.append(probe)
        return ans
    
    def print_fsa_from_idfile(self,filename):
        ids = self.ids_from_file(filename)
        self.print_fsa_from_ids(ids)

    def seqs_from_idfile(self,filename,QUIET=''):
        return(self.seqs_from_ids(self.ids_from_file(filename,QUIET)))

    def seqs_from_ids(self,ids):
        seqs = []
        f_ids = self.filter(ids)
        for id in f_ids:
            seqs.append(self.probes[id])
        return(seqs)

    def reduce_probelist(self,ids):
        """
        ProbeSet.reduce_probelist(ids) --
                 Change the the internal probelist.  Keep only sequences associated with the supplied IDs.
        """
        newprobes = {}
        for id in ids:
            newprobes[id] = self.probes[id]
        self.probes = newprobes

def AmbigToRegExp(expression):
    """
    AmbigToRegExp(expression)
             Converts a string IUPAC ambiguity codes into a regular expressions. 
    """
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

def permute(letters, depth, seqs=[''],curdepth=0):
    """
    permute(letters, depth, seqs=[''],curdepth=0) --
           Generate all possible sequences from the alphabet
           "letters" of length "depth" (sorry about the
           variable names).  For example, "permute('ACGT',4)"
           generates all 256 possible 4-letter DNA sequences.

           Duplicated in TAMO.util.PermuteTools
           
    """
    newseqs = []
    for seq in seqs:
        for letter in letters:
            newseqs.append(seq + letter)
    if depth > curdepth:
        return(permute(letters,depth,newseqs,curdepth + 1))
    else:
        return(seqs)

	
revcomp_memo  = {}
revcomp_trans = string.maketrans("AGCTagctnNRYWSMKHBVD", "TCGAtcganNYRWSKMDVBH")
def revcomp(dna):
    """
    revcomp(dna) --
           return reverse complement of a DNA sequence
    """
    try: ans = revcomp_memo[dna]
    except:
        comp = dna.translate(revcomp_trans)
        lcomp = list(comp)
        lcomp.reverse()
        ans = ''.join(lcomp)
        revcomp_memo[dna] = ans
    return ans
    
def top_nmers_fsa(N,filename):
    """
    top_nmers_fsa(N,filename) --
            Return the n-mers (or k-mers) of length N in sequencs in the fasta-formatted
            input file.
    """
    seqs = Fasta.seqs(filename)
    return(top_nmers_seqs(N,seqs))
            
def top_nmers_seqs(N,seqs):
    """
    top_nmers_seqs(N,seqs) --
            Return the n-mers (or k-mers) that occur at in at least approximately 10% of the input
            sequences.
    """
    nmers = []
    mincount = max(2, 0.1*float(len(seqs)))
    for nmer,count in MotifTools.top_nmers(N,seqs,'tuples'):
        numN = N - len(nmer.replace('N',''))
        if count > mincount  and  numN+5 <= len(nmer):
            nmers.append(nmer)
    return(nmers)

def gapped_motifs(seq,maxgap=12):
    """
    gapped_motifs(seq,maxgap=12) -- 
         Takes input text string, say CGG, and builds a list of all tandem, inverted, and everted repeats
         with "N" gaps ranging from 0 to maxgap (12). [CGGCCG CGGCGG CCGCGG, CGGnCCG CGGnCGG CCGnCGG, ....]
    """
    seqrc = revcomp(seq)
    motifs = []
    for l,r in [(seq,seq), (seq, seqrc), (seqrc, seq)]:
        for i in range(maxgap):
            motifs.append('%s%s%s'%(l, 'N'*i, r))
    return(motifs)


if __name__ == '__main__': main()
