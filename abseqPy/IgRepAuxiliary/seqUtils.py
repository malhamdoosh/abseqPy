'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from __future__ import print_function

import matplotlib
matplotlib.use('agg')
import gc
import sys
import pickle
import random
import bisect
import matplotlib.pyplot as plt
import os

from os.path import exists
from collections import Sequence, defaultdict
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio import SeqIO, motifs, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet
from Bio import Alphabet

from abseqPy.IgRepertoire.igRepUtils import alignListOfSeqs, safeOpen, detectFileFormat
from abseqPy.config import WEBLOGO
from abseqPy.logger import printto, LEVEL
from abseqPy.utilities import ShortOpts, requires, quote


# the following are conditionally imported in functions that require them to reduce abseq's dependency list
# It's here for a simple glance of required dependencies
# from TAMO.Clustering.UPGMA import UPGMA
# from TAMO.Clustering.UPGMA import DFUNC
# from TAMO.Clustering.UPGMA import print_tree
# from TAMO.Clustering.UPGMA import create_tree_phylip
# from TAMO.Clustering.UPGMA import print_tree_id
# from TAMO import MotifTools
# from TAMO.MotifTools import Motif


def readSeqFileIntoDict(seqFile, outDict=None, stream=None):
    printto(stream, "Processing {} ... loading sequences into dictionary".format(os.path.basename(seqFile)))
    format = detectFileFormat(seqFile)
    if outDict is None:
        outDict = {}
    with safeOpen(seqFile) as fp:
        for rec in SeqIO.parse(fp, format):
            outDict[rec.id] = str(rec.seq)
    return outDict


def generateMotif(sequences, name, alphabet, filename, 
                  align=False, transSeq=False, protein=False, weights=None, outDir=None,
                  threads=2, stream=None):
    """

    :param sequences: list of strings
                    list of sequences used to find motifs

    :param name: string
                    sample name

    :param alphabet: Bio.Alphabet

    :param filename: string
                    output filename

    :param align: bool
                    use CLUSTAL OMEGA for sequence alignment?

    :param transSeq: bool
                    should sequences be translated to protein?

    :param protein: bool

    :param weights: list
                    list of weights for provided sequences

    :param outDir: string
                    path to output directory

    :param threads: int
                    number of threads to use

    :param stream: logger
                    calls printto stream object

    :return: BioPy.Motif object
    """
    if exists(filename):
        printto(stream, "\t" + name + " motif logo was found", LEVEL.WARN)
        return

    # check whether sequences should be translated                 
    if transSeq:
        seqs = []               
        for rec in sequences:
            seq = Seq(rec[: len(rec) - (len(rec) % 3)]).translate(to_stop=False)
            seqs.append(str(seq))
    else:
        seqs = sequences

    # sample sequences if they are too many  
    if (len(seqs) > 1*10**5 or
            (weights is not None and sum(weights) > 1*10**5)):
        random.seed(1986)
        seqs = weightedSample(seqs, weights, int(1*10**5))
        if align:
            seqs = random.sample(seqs, 10000)

    # perform multiple sequence alignment on a sample of 10000 sequences 
    if align and len(seqs) > 1:
        alignedSeq = alignListOfSeqs(seqs, outDir, threads=threads, name=name, stream=stream)
    else:
        # if alignment is not required, add "-" to short sequences
        L = map(len, seqs)
        if min(L) != max(L):
            # print('\t\t- is being added to short sequences ...[%d, %d[' % (min(L), max(L)))
            if '-' not in alphabet.letters:
                alphabet.letters += '-'
            alignedSeq = []
            m = max(L)
            for s in seqs:
                if len(s) < m:
                    alignedSeq.append(s + '-'*(m-len(s)))
        else:
            alignedSeq = seqs

    # create the sequence motif and encode it into PFM
    printto(stream, "\tMotif logo is being created for %s ..." % name)
    m = motifs.create(alignedSeq, alphabet)  # print(m.counts)
    # create sequence logo
    generateMotifLogo(m, filename, outDir, not transSeq and not protein, stream=stream)
    return m
    
    
def createAlphabet(align=False, transSeq=False, extendAlphabet=False, protein=False):
    if not transSeq and not protein:
        alphabet = Alphabet.DNAAlphabet()
        alphabet.letters = "ACGT" if not extendAlphabet else "ACGTN"
    else:
        alphabet = Alphabet.ProteinAlphabet()
        alphabet.letters = IUPACProtein().letters + ('*' if not extendAlphabet else '*X')
    if align:
        alphabet.letters += '-'
    return alphabet


@requires("TAMO")
def generateMotifs(seqGroups, align, outputPrefix, transSeq=False,
                        extendAlphabet=False, clusterMotifs=False, protein=False, threads=2, stream=None):
    from TAMO.MotifTools import Motif
    ighvMotifs = []
    if clusterMotifs and 'gene' in outputPrefix:
        findMotifClusters(ighvMotifs, outputPrefix, stream=stream)
    printto(stream, '\t\tPWMs, consensus and logos are being generated for {} motifs ... '.format(len(seqGroups)))
    pwmFile = open(outputPrefix + '_pwm.txt', 'w')
    consensusFile = open(outputPrefix + '_consensus.txt', 'w')
    logosFolder = outputPrefix + '_logos'

    if not os.path.exists(logosFolder):
        os.makedirs(logosFolder)

    # create the sequence alphabet: DNA or Protein
    alphabet = createAlphabet(align, transSeq, extendAlphabet, protein)
    groups = seqGroups.keys()
    groups.sort()        
    
    for group in groups:    
        filename = os.path.join(logosFolder, group.replace('/', '') + '.png')
        seqs = seqGroups[group]
        m = generateMotif(seqs, group, alphabet, filename, align, transSeq, protein, outDir=logosFolder,
                          threads=threads, stream=stream)
        if m is None:
            # motif file found, no further work required
            return
        motifSeqs = m.instances
        pwm = m.counts.normalize(pseudocounts=None)  # {'A':0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}
        consensusMax = str(m.consensus)      
               
        pwmFile.write('#{} {} sequences\n'.format(group, len(motifSeqs)))
        pwmFile.write(str(pwm))  
        consensusFile.write('>{} max_count\n'.format(group))
        consensusFile.write(consensusMax + '\n')      
    #             print(str(m.anticonsensus)) # smallest values in the columns
        if not transSeq and not align and not protein:
            consensusIupac = str(m.degenerate_consensus)
    #             print(consensusIupac) # IUPAC ambiguous nucleotides            
            consensusFile.write('>{} degenerate\n'.format(group))
            consensusFile.write(consensusIupac + '\n')
        
        pwmFile.flush()
        consensusFile.flush()
        gc.collect()
        if clusterMotifs and len(motifSeqs) > 10:
            motif = Motif(map(lambda x: str(x), motifSeqs),
                          backgroundD={'A': 0.6, 'C': 0.4, 'G': 0.4, 'T': 0.6}, id=group)
            motif.addpseudocounts(0.1)
            ighvMotifs.append(motif)
            
    pwmFile.close()
    consensusFile.close()      
    gc.collect()
    printto(stream, "\tPosition weight matrices are written to " + os.path.basename(outputPrefix + '_pwm.txt'))
    printto(stream, "\tConsensus sequences are written to " + os.path.basename(outputPrefix + '_consensus.txt'))
    if clusterMotifs:
        findMotifClusters(ighvMotifs, outputPrefix, stream=stream)


@requires("TAMO")
def findMotifClusters(ighvMotifs, outputPrefix, stream=None):
    from TAMO.Clustering.UPGMA import UPGMA
    from TAMO.Clustering.UPGMA import DFUNC
    from TAMO.Clustering.UPGMA import print_tree_id
    # cluster using a variant of the UPGMA algorithm implemented in the TAMO package
    
    motifsFile = os.path.abspath(outputPrefix + '_motifs.tamo')
    if not exists(motifsFile):
        if len(ighvMotifs) > 0:
            pickle.dump(ighvMotifs, open(motifsFile, 'wb'))            
    else:
        ighvMotifs = pickle.load(open(motifsFile, 'rb'))

    prefixName, sampleName = os.path.split(outputPrefix)
    dendrogramDirectory = os.path.join(prefixName, 'dendrograms')
    if not exists(dendrogramDirectory):
        os.makedirs(dendrogramDirectory)

    if len(ighvMotifs) > 0:
        groupedMotifs = defaultdict(list)
        for m in ighvMotifs:
            ighv = m.id.split('-')[0].split('/')[0]
            groupedMotifs[ighv].append(m)
        try:
            motifClustersFile = os.path.join(dendrogramDirectory, sampleName + '_pwm_clusters.txt')

            _old_stdout = sys.stdout
            sys.stdout = open(motifClustersFile, 'w')

            for ighv in groupedMotifs.keys():
                newickdendrogramFile = os.path.join(dendrogramDirectory, sampleName + '_{}_newick.dnd'.format(ighv))
                tree = UPGMA(groupedMotifs[ighv], DFUNC)
                print_tree_id(tree)

                saveNewickdendrogram(newickdendrogramFile, tree, sys.stdout, title=(ighv + " family clustering"), logger=stream)

            lists = groupedMotifs.values()
            tree = UPGMA([m for lst in lists for m in lst], DFUNC)
            print_tree_id(tree)

            newickdendrogramFile = os.path.join(dendrogramDirectory, sampleName + '_newick.dnd')
            saveNewickdendrogram(newickdendrogramFile, tree, sys.stdout, title="Clustering of all IGHV", logger=stream)

            sys.stdout.close()
            sys.stdout = _old_stdout

            printto(stream, "\tMotif clusters were written to " + os.path.basename(motifClustersFile))
        except Exception as e:
            printto(stream, "Motifs couldn't be clustered! Error: {}".format(str(e)), LEVEL.ERR)


def saveNewickdendrogram(newickClusterFile, tree, stream, title="", logger=None):
    """
    :param newickClusterFile:
    :param tree:  UPGMA object
    :param stream:
    :param title:
    :param logger:
    :return:
    """
    from TAMO.Clustering.UPGMA import create_tree_phylip
    desc = '' if not title else " for {} ".format(title)

    # get phylip newick syntax
    phylipTree = create_tree_phylip(tree)
    with open(newickClusterFile, 'w') as newickfp:
        newickfp.write(phylipTree)

    printto(logger, "Newick dendrogram{}written to ".format(desc) + os.path.basename(newickClusterFile))

    # show ascii art
    phylipTree = Phylo.read(newickClusterFile, format='newick')

    try:
        print("\n\nASCII phylip tree{}:\n".format(desc), file=stream)
        Phylo.draw_ascii(phylipTree, file=stream)
    except ZeroDivisionError:
        # if the weights are 0
        print("\t Not drawn because of 0 weights", file=stream)
        pass

    # plot dendrogram in matplotlib
    phylipTree.ladderize()
    fig, axes = plt.subplots(figsize=(8, 5))
    Phylo.draw(phylipTree, do_show=False, axes=axes, show_confidence=True)
    axes.set_title(title)
    fig.savefig(newickClusterFile.replace('.dnd', '.png'), dpi=300)
    plt.close()


@requires('weblogolib')
def generateMotifLogo(m, filename, outdir='.', dna=True, stream=None):
    instances = m.instances
    records = []

    for i in range(len(instances)):
        records.append(SeqRecord(instances[i], id=str(i)))

    # guarantees tha tmp file is unique (in the case where multiprocessing might override the file)
    basename = os.path.splitext(os.path.basename(filename))[0]
    tmpFile = os.path.join(outdir, 'temp_seq_logos_{}.fasta'.format(basename))

    SeqIO.write(records, tmpFile, 'fasta')

    weblogo = ShortOpts(exe=WEBLOGO, f=quote(tmpFile), o=quote(filename), A=("dna" if dna else "protein"), F="png",
                        n=200, D="fasta", s="medium", c=("classic" if dna else "auto"),
                        X="NO", Y="NO")\
        .append("--errorbars NO --fineprint CSL --resolution 600")
    # printto(stream, "Executing " + str(weblogo))
    weblogo()
    os.remove(tmpFile)
        

def maxlen(x):
    return max(map(len, x))


def weightedSample(population, weights, k):
    if weights is not None:
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

