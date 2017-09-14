'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

import  IgRepertoire
from Bio.Seq import  Seq
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio.Alphabet import Alphabet
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
from collections import Sequence
import bisect
from numpy import dtype
#from IgRepertoire.igRepUtils import alignListOfSeqs


def readSeqFileIntoDict(seqFile, format = "fastq", outDict = None):
    if (outDict is None):
        outDict = {}
    try:
        if format == "fastq":
            with open(seqFile) as h:
                while True:
                    line = h.readline()
                    if not line:
                        break
                    id = line.strip("\n")[1:].split()[0]
                    seq = h.readline().strip("\n") 
                    outDict[id] = seq
                    h.readline()
                    h.readline()
#                     print(id)
        elif format == "fasta":
            with open(seqFile) as h:
                while True:            
                    line = h.readline()
                    if not line:
                        break
                    line = line.strip("\n")
                    if line.startswith(">"):
                        id = line[1:].split()[0]
                        outDict[id] = ""
                    else:
                        outDict[id] += line                            
        else:
            raise Exception("Unknown sequence file format")
    except Exception as e: 
        print("Something went wrong while reading a sequence file")
        raise e
    return outDict

def generateMotif(sequences, name, alphabet, filename, 
                  align = False, transSeq = False, protein =False, weights = None, outDir=None):
    if (exists(filename)):
        print("\t" + name + " motif logo was found" )
        return        
    # check whether sequences should be translated                 
    if transSeq:
        seqs = []               
        for rec in sequences:
            seq = Seq(rec).translate(to_stop=False)                   
            seqs.append(str(seq))
    else:
        seqs = sequences
    # sample sequences if they are too many  
    if (len(seqs) > 1*10**5 or 
        ( weights is not None and sum(weights) > 1*10**5)):
        random.seed(1986)
#         print(sum(weights))
        seqs = weightedSample(seqs, weights, int(1*10**5)) 
        if align:
            seqs = random.sample(seqs, 10000) 
    # perform multiple sequence alignment on a sample of 10000 sequences 
    if align and len(seqs) > 1:
        alignedSeq = IgRepertoire.igRepUtils.alignListOfSeqs(seqs, outDir)
#                 print(alignedSeq[:10])
    else:                
        # if alignment is not required, add "-" to short sequences
        L = map(len, seqs)
        if (min(L) != max(L)):
            #print('\t\t- is being added to short sequences ...[%d, %d[' % (min(L), max(L)))
            if '-' not in alphabet.letters: alphabet.letters += '-'
            alignedSeq = []                    
            m = max(L)
            for s in seqs:
                if len(s) < m:
                    alignedSeq.append(s + '-'*(m-len(s)))
        else:
            alignedSeq = seqs    
    # create the sequence motif and encode it into PFM
    print("\tMotif logo is being created for %s ..." % (name))
    m = motifs.create(alignedSeq, alphabet)  #             print(m.counts)
    # create sequence logo
    generateMotifLogo(m, filename, outDir, not transSeq and not protein)
    return m
    
    
def createAlphabet(align = False, transSeq=False, 
                   extendAlphabet = False, protein = False):
    if not transSeq and not protein:
        alphabet = Alphabet()
        alphabet.letters = "ACGT" if not extendAlphabet else "ACGTN"            
    else:
        alphabet = IUPACProtein()
        alphabet.letters += '*' if not extendAlphabet else '*X'         
    if align:
        alphabet.letters += '-'  
    return alphabet

def generateMotifs(seqGroups, align, outputPrefix, transSeq=False,
                        extendAlphabet=False, clusterMotifs=False, protein=False):  
    ighvMotifs = []
    if (clusterMotifs and 'gene' in outputPrefix):
        findMotifClusters(ighvMotifs, outputPrefix)                
    print("\t\tPWMs, consensus and logos are being generated for %d motifs ... " % (len(seqGroups)))      
    pwmFile = open(outputPrefix + '_pwm.txt', 'w')
    consensusFile = open(outputPrefix + '_consensus.txt', 'w')
    logosFolder = outputPrefix + '_logos/'
    os.system('mkdir ' + logosFolder)
    # create the sequence alphabet: DNA or Protein
    alphabet = createAlphabet(align, transSeq, extendAlphabet, protein)
    groups = seqGroups.keys()
    groups.sort()        
    
    for group in groups:    
        filename = logosFolder + group.replace('/', '') + '.png'    
        seqs = seqGroups[group]
        m = generateMotif(seqs, group, alphabet, filename, align, transSeq, 
                          protein, outDir=logosFolder)
        motifSeqs = m.instances
        pwm = m.counts.normalize(pseudocounts=None)  # {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
        consensusMax = str(m.consensus)      
               
        pwmFile.write('#%s %d sequences\n' % 
                      (group, len(motifSeqs)))
        pwmFile.write(str(pwm))  
        consensusFile.write('>%s max_count\n' % (group))
        consensusFile.write(consensusMax + '\n')      
    #             print(str(m.anticonsensus)) # smallest values in the columns
        if (not transSeq and not align and not protein):
            consensusIupac = str(m.degenerate_consensus)
    #             print(consensusIupac) # IUPAC ambiguous nucleotides            
            consensusFile.write('>%s degenerate\n' % (group))
            consensusFile.write(consensusIupac + '\n')
        
        pwmFile.flush()
        consensusFile.flush()
        gc.collect()
        if (clusterMotifs and len(motifSeqs) > 10):
            motif = Motif(map(lambda x: str(x), motifSeqs), 
                     backgroundD={'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6},
                     id = group)
            motif.addpseudocounts(0.1)
            ighvMotifs.append(motif)
            
    pwmFile.close()
    consensusFile.close()      
    gc.collect()
    print("\tPosition weight matrices are written to " + outputPrefix + '_pwm.txt')
    print("\tConsensus sequences are written to " + outputPrefix + '_consensus.txt')
    if (clusterMotifs):
        findMotifClusters(ighvMotifs, outputPrefix)
        
        
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

def generateMotifLogo(m, filename, outdir, dna=True):
    instances = m.instances
    records = []
    for i in range(len(instances)):
        records.append(SeqRecord(instances[i], id=`i`))
    tmpFile = (((outdir+'/') if outdir else "") + 'temp_seq_logos.fasta').replace('//', '/')
    SeqIO.write(records, tmpFile, 'fasta')

    command = "%s -f %s  -o %s -A %s -F png -n 200 -D fasta -s medium "
    command += "-c %s --errorbars NO --fineprint CSL --resolution 600 -X NO -Y NO"

    os.system(command % (WEBLOGO, tmpFile,
                         filename, "dna" if dna else "protein",
                         "classic" if dna else "auto"))
    os.remove(tmpFile)
        

def maxlen(x):
    return max(map(len, x))

#TODO: Look for faster and ACCURATE weighted sampling approach 
def weightedSampleFast(population, weights, k):
    if (weights is not None):
        from fast_sampler import FastSampler
        from numpy import array
        weights = array(weights, dtype='d')        
        h = FastSampler(len(population), max(weights), min(weights))
        for i in range(len(population)):
            h.add(i, weights[i])
        s = map(lambda x : population[h.sample()], range(k))        
        return s
    else:
        return random.sample(population, k)
    
def weightedSample(population, weights, k):
    if (weights is not None):
        return random.sample(WeightedPopulation(population, weights), k)
    else:
        return random.sample(population, k)
# from http://stackoverflow.com/questions/13047806/weighted-random-sample-in-python
class WeightedPopulation(Sequence):
    def __init__(self, population, weights):
        assert len(population) == len(weights) > 0
        self.population = population
        self.cumweights = []
        cumsum = 0 # compute cumulative weight
        for w in weights:
            cumsum += w   
            self.cumweights.append(cumsum)  
    def __len__(self):
        return self.cumweights[-1]
    def __getitem__(self, i):
        if not 0 <= i < len(self):
            raise IndexError(i)
        return self.population[bisect.bisect(self.cumweights, i)]
    
    
    
