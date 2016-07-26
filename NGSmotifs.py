from Bio.Seq import  Seq
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio.Alphabet import Alphabet
from igRepUtils import alignListOfSeqs
from os.path import exists
import os
from Bio import SeqIO, motifs
from Bio.SeqRecord  import SeqRecord
from config import WEBLOGO
from TAMO.Clustering.UPGMA import UPGMA
from TAMO.Clustering.UPGMA import DFUNC
from TAMO.Clustering.UPGMA import print_tree
from TAMO.MotifTools import Motif
import gc
from TAMO.Clustering.UPGMA import print_tree_id
from TAMO import MotifTools
import sys
import pickle
from TAMO.Clustering.UPGMA import create_tree_phylip
import random
import numpy
from matplotlib.pyplot import *
from collections import Counter



def findMotifClusters(ighvMotifs, outputPrefix):
    # cluster using a variant of the UPGMA algorithm implemented in the TAMO package 
    
    motifsFile = os.path.abspath(outputPrefix + '_motifs.tamo')
    if (not exists(motifsFile)):
        if (len(ighvMotifs) > 0):
            pickle.dump(ighvMotifs, open(motifsFile, 'wb'))            
    else:
        ighvMotifs = pickle.load(open(motifsFile, 'rb'))
#         print(ighvMotifs)
    if (len(ighvMotifs) > 0): 
        groupedMotifs = {}
        for m in ighvMotifs:
            ighv = m.id.split('-')[0].split('/')[0]
            if (groupedMotifs.get(ighv, None) is None):
                groupedMotifs[ighv] = []
            groupedMotifs[ighv].append(m)
        try:
            motifClustersFile = outputPrefix + '_pwm_clusters.txt'            
            _old_stdout = sys.stdout
            sys.stdout = open(motifClustersFile, 'w')
            for ighv in groupedMotifs.keys():
                tree = UPGMA(groupedMotifs[ighv], DFUNC)
#                 print_tree(tree)
                print_tree_id(tree)
                print(create_tree_phylip(tree))
                sys.stdout.flush()
            lists = groupedMotifs.values()
            tree = UPGMA([m for list in lists for m in list], DFUNC)
#                 print_tree(tree)
            print_tree_id(tree)
            print(create_tree_phylip(tree))
            sys.stdout.close()
            sys.stdout = _old_stdout
            print("\tMotif clusters were written to " + motifClustersFile)
        except:
            print("Motifs couldn't be clustered!")
#             raise


def generateMotifs(ighvSignals, align, outputPrefix, transSeq=False,
                        extendAlphabet=False, clusterMotifs=False, protein=False):  
        ighvMotifs = []
        if (clusterMotifs and 'gene' in outputPrefix):
            findMotifClusters(ighvMotifs, outputPrefix)                
        print("\t\tPWMs, consensus and logos are being generated for %d motifs ... " % (len(ighvSignals)))      
        pwmFile = open(outputPrefix + '_pwm.txt', 'w')
        consensusFile = open(outputPrefix + '_consensus.txt', 'w')
        logosFolder = outputPrefix + '_logos/'
        os.system('mkdir ' + logosFolder)
        # create the sequence alphabet: DNA or Protein
        if not transSeq and not protein:
            alphabet = Alphabet()
            alphabet.letters = "ACGT" if not extendAlphabet else "ACGTN"            
        else:
            alphabet = IUPACProtein()
            alphabet.letters += '*' if not extendAlphabet else '*X'         
        if align:
            alphabet.letters += '-'  
        ighvs = ighvSignals.keys()
        ighvs.sort()        
        
        for ighv in ighvs:    
            filename = logosFolder + ighv.replace('/', '') + '.png'    
            if (exists(filename)):
                print("File found ... " + filename.split("/")[-1])
                continue        
            # check whether sequences should be translted                 
            if transSeq:
                signals = []               
                for rec in ighvSignals[ighv]:
                    seq = Seq(rec).translate(to_stop=False)                   
                    signals.append(str(seq))
            else:
                signals = ighvSignals[ighv]  
            if (len(signals) > 1*10**5):
                random.seed(1986)
                signals = random.sample(signals, int(1*10**5)) 
                if align:
                    signals = random.sample(signals, 10000) 
            # perform multiple sequence alignment on a sample of 5000 sequences 
            if align and len(signals) > 1:
                alignedSeq = alignListOfSeqs(signals)
#                 print(alignedSeq[:10])
            else:                
                # if alignment is not required, add - to short sequences
                L = map(len, signals)
                if (min(L) != max(L)):
                    print('\t\t\t- is being added to short sequences ...[%d, %d[' % (min(L), max(L)))
                    if '-' not in alphabet.letters: alphabet.letters += '-'
                    alignedSeq = []                    
                    m = max(L)
                    for s in signals:
                        if len(s) < m:
                            alignedSeq.append(s + '-'*(m-len(s)))
                else:
                    alignedSeq = signals    
            # create the sequence motif and encode it into PFM
            print("\t\t\tMotif is being created for %s ..." % (ighv))
            m = motifs.create(alignedSeq, alphabet)  #             print(m.counts)
            pwm = m.counts.normalize(pseudocounts=None)  # {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
            if (clusterMotifs and len(alignedSeq) > 10):
                ighvMotifs.append(Motif(alignedSeq, 
                                         backgroundD={'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6},
                                         id = ighv))
                ighvMotifs[-1].addpseudocounts(0.1)
#             if (len(ighvMotifs) > 10):
#                 break
#                 print(ighvMotifs[ighv])            
            pwmFile.write('#%s %d sequences\n' % 
                          (ighv, len(alignedSeq)))
            pwmFile.write(str(pwm))  
#             print(pwm)         
            consensusMax = str(m.consensus)        
#             print(consensusMax) # largest values in the columns
            consensusFile.write('>%s max_count\n' % (ighv))
            consensusFile.write(consensusMax + '\n')      
#             print(str(m.anticonsensus)) # smallest values in the columns
            if (not transSeq and not align and not protein):
                consensusIupac = str(m.degenerate_consensus)
    #             print(consensusIupac) # IUPAC ambiguous nucleotides            
                consensusFile.write('>%s degenerate\n' % (ighv))
                consensusFile.write(consensusIupac + '\n')
            # create sequence logo
            generateMotifLogo(m, filename,
                              not transSeq and not protein)
            pwmFile.flush()
            consensusFile.flush()
            gc.collect()
#             m.weblogo(logosFolder + ighv.replace('/', '') + '.png')
#             if (consensusMax != consensusIupac):
#                 print(consensusIupac, consensusMax)
                
        pwmFile.close()
        consensusFile.close()      
        gc.collect()
        print("\tPosition weight matrices are written to " + outputPrefix + '_pwm.txt')
        print("\tConsensus sequences are written to " + outputPrefix + '_consensus.txt')
        if (clusterMotifs):
            findMotifClusters(ighvMotifs, outputPrefix)
        
        
def generateMotifLogo(m, filename, dna=True):
    instances = m.instances
    records = []
    for i in range(len(instances)):
        records.append(SeqRecord(instances[i], id=`i`))
    SeqIO.write(records, 'temp_seq_logos.fasta', 'fasta')
    
    command = "%s -f %s  -o %s -A %s -F png -n 200 -D fasta -s medium " 
    command += "-c %s --errorbars NO --fineprint CSL --resolution 600 -X NO -Y NO" 
    os.system(command % (WEBLOGO, 'temp_seq_logos.fasta', 
                         filename, "dna" if dna else "protein",
                         "classic" if dna else "auto"))

    os.system('rm temp_seq_logos.fasta')     
        
      
AA = ["GAST", "CVILPFYMW", "NQH", "DE", "KR"]
AA_colours = numpy.concatenate((
    cm.Oranges( (1+numpy.arange(len(AA[0]), dtype=float)) / (len(AA[0])+1)),
    cm.Greens(  (1+numpy.arange(len(AA[1]), dtype=float)) / (len(AA[1])+1)),
    cm.Purples( (1+numpy.arange(len(AA[2]), dtype=float)) / (len(AA[2])+1)),
    cm.Reds(    (1+numpy.arange(len(AA[3]), dtype=float)) / (len(AA[3])+1)),
    cm.Blues(   (1+numpy.arange(len(AA[4]), dtype=float)) / (len(AA[4])+1))
))

AA = ''.join(AA)  
        
def barLogo(counts, _title, filename):
    if (exists(filename)):
        print("File found ... " + filename.split('/')[-1])
        return
    fig, ax = subplots(figsize=(8,5))
    bar_fractions = [ [ ct.get(aa, 0) / float(sum(ct.values())) for aa in AA ] for ct in counts ]
    by_aa      = [ [] for aa in AA ]
    by_aa_base = [ [] for aa in AA ]
    for bf in bar_fractions:
        s = 0.0
        for i,aa in enumerate(AA):
            by_aa[i].append(bf[i])
            by_aa_base[i].append(s)
            s += bf[i]

    ax.set_title(_title, fontsize=20)
    for i,aa in enumerate(AA):
        ax.bar(numpy.arange(len(bar_fractions)) + .05, by_aa[i], 
               width=0.9, bottom=by_aa_base[i], color=AA_colours[i], 
               label=AA[i], lw=0)   
    ax.set_ylim(0,1)
    ax.set_xticks(numpy.arange(len(counts))+.5 )
    ax.set_xticklabels([ ct.most_common(1)[0][0] for ct in counts ])
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.legend(loc = 'upper right', bbox_to_anchor = (1.1, 1),fontsize='x-small')
    fig.savefig(filename, dpi=300)
    
def maxlen(x):
    return max(map(len, x))

def generateProteinLogos(proteinSeqs, outputPrefix):
    logosFolder = outputPrefix + '_logos/'    
    groups = proteinSeqs.keys()
    groups.sort()        
    if (exists(logosFolder + groups[-1] + '.png')):
        print("Logos were already generated.")
        return
    print("\tProtein sequence logos are being generated .... ")   
    os.system('mkdir ' + logosFolder) 
    for group in groups: 
        print("\t\t" + group)
        seqs = proteinSeqs[group]
        m = maxlen(seqs)
        if m > 32:
            m = 32
        aa_counts = [ Counter(c[x] for c in seqs if len(c) > x) for x in range(m)]
#         print(aa_counts)
        barLogo(aa_counts, "{} ({:,})".format(group.upper(), len(seqs)), logosFolder + group + '.png')
    
    
    
    
    