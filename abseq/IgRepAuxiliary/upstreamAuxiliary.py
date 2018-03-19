'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
import os
import gc

from numpy import Inf
from Bio import SeqIO

from abseq.IgRepReporting.igRepPlots import plotSeqLenDist, plotSeqLenDistClasses, plotDist
from abseq.IgRepertoire.igRepUtils import gunzip, writeCountsToFile, compressCountsFamilyLevel, compressCountsGeneLevel
from abseq.logger import LEVEL, printto


def plotUpstreamLenDist(upstreamFile, expectLength, name, stream=None):
    """
    plots length distribution of upstream sequences. 4 different sets of plots are generated:
        1. class-level distribution
        2. gene-level distribution
        IF expectLength != INF:
            3. class-level distribution for upstream seqs that are shorter than expectedLength
            4. gene-level distribution for upstream seqs that are shorter than expectedLength

    :param upstreamFile:
                upstream sequences (FASTA file)

    :param expectLength:
                expected length of upstream sequences (int)

    :param name:
                sample name (string)

    :param stream:
                logging stream

    :return:
                None
    """
    if not upstreamFile.endswith(".fasta"):
        raise Exception("Expected {} to be a FASTA file".format(os.path.basename(upstreamFile.rstrip(os.path.sep))))

    fileExt = 'fasta'

    outputFile = upstreamFile.replace('.fasta', '_dist.png')
    plotSeqLenDist(upstreamFile, name, outputFile, fileExt, stream=stream)

    outputFile = upstreamFile.replace('.fasta', '_dist_class.png')
    plotSeqLenDistClasses(upstreamFile, name, outputFile, fileExt, stream=stream)

    # if user provided values to upstream (and it's not Inf)
    if expectLength != Inf:
        # plot seqs that do not meet the expectedLength (shorter than expected length)
        outputFile = upstreamFile.replace('.fasta', '_dist_short.png')
        plotSeqLenDist(upstreamFile, name, outputFile, fileExt, expectLength - 1, stream=stream)

        outputFile = upstreamFile.replace('.fasta', '_dist_short_class.png')
        plotSeqLenDistClasses(upstreamFile, name, outputFile, fileExt, expectLength - 1, stream=stream)


def extractUpstreamSeqs(cloneAnnot, recordFile, upstream, upstreamFile, stream=None):
    """
    extract the upstream DNA sequences and write them into a FASTA file named upstreamFile
    :param cloneAnnot:
                cloneAnnot DataFrame

    :param recordFile:
                raw record file (string)

    :param upstream:
                list of 2 numbers, denoting [start, end] inclusive in 1-index. np.Inf is also allowed for end value

    :param upstreamFile:
                output FASTA filename

    :param stream:
                logging stream object

    :return:
                None
    """
    printto(stream, "\tExtracting the upstream sequences ... ")

    # alignments with - strand
    revAlign = 0
    # num. seqs with trimmed beginning (vstart > 3)
    trimmedBegin = 0
    # num. seqs with sequences shorter than expected upstream length (len(seq) < expectLength)
    # @see expectLength
    trimmedUpstream = 0
    # excluded sequences because end <= 1
    noSeq = 0
    # num. processed sequences
    procSeqs = 0
    # buffer to hold sequences before flushing into file
    recordsBuffer = []
    # max buffer size allowed
    maxBufferSize = int(10.0 ** 5) / 2

    # expected upstream length = expectLegth (end - start + 1) where start,end are both 1-index
    expectLength = upstream[1] - upstream[0] + 1
    queryIds = cloneAnnot.index

    # NOTE: SeqIO.index can only index string filenames and it has to be unzipped
    _, ext = os.path.splitext(os.path.basename(recordFile.rstrip(os.path.sep)))
    records = SeqIO.index(gunzip(recordFile), ext.lstrip('.'))

    with open(upstreamFile, 'w') as fp:
        for id_ in queryIds:
            record = records[id_]
            qsRec = cloneAnnot.loc[record.id]
            if qsRec.strand != 'forward':
                revAlign += 1
                record.seq = record.seq.reverse_complement()
            if qsRec.vstart <= 3:
                end = qsRec.vqstart - upstream[0] - qsRec.vstart + 1
                if end <= 1:
                    noSeq += 1
                else:
                    start = max(1, qsRec.vqstart - upstream[1] - qsRec.vstart + 1)
                    record.seq = record.seq[int(start - 1):int(end)]
                    if expectLength != Inf and len(record.seq) < expectLength:
                        trimmedUpstream += 1
                    record.id = record.id + '|' + qsRec.vgene
                    record.description = ""
                    recordsBuffer.append(record)
                    procSeqs += 1
                    if procSeqs % maxBufferSize == 0:
                        printto(stream, '{}/{} sequences have been processed ... '.format(procSeqs, len(queryIds)))
                        SeqIO.write(recordsBuffer, fp, 'fasta')
                        recordsBuffer = []
            else:
                trimmedBegin += 1

        # flush remaining sequences
        if len(recordsBuffer) > 0:
            printto(stream, '{}/{} sequences have been processed ... '.format(procSeqs, len(queryIds)))
            SeqIO.write(recordsBuffer, fp, 'fasta')

    if revAlign > 0:
        printto(stream,
                "\t\t\t{} sequences are in reversed alignment ... ".format(revAlign),
                LEVEL.INFO)

    if trimmedBegin > 0:
        printto(stream, "\t\t\tThe query sequence is not aligned within 3bp of the IGV start "
                        "position ... {} found and excluded!".format(trimmedBegin), LEVEL.WARN)

    if trimmedUpstream > 0:
        printto(stream,
                "\t\t\tUpstream sequences shorter than the expected length are detected ... {} found"
                .format(trimmedUpstream),
                LEVEL.WARN)

    if noSeq > 0:
        printto(stream,
                "\t\t\tNo upstream sequence can be extracted (too short) for {} sequences.".format(noSeq),
                LEVEL.WARN)
    gc.collect()


# todo Mon Mar 19 18:13:58 AEDT 2018
def loadValidSequences():
    pass


# todo Mon Mar 19 18:13:58 AEDT 2018
def analyzeSequences():
    pass

def writeCountsCategoriesToFile(countsVariant, sampleName, filePrefix, title=''):
    writeCountsToFile(countsVariant,
                      filePrefix + 'variant.csv')
    # gene level
    countsVariant = compressCountsGeneLevel(countsVariant)
    writeCountsToFile(countsVariant,
                      filePrefix + 'gene.csv')
    plotDist(countsVariant, sampleName,
             filePrefix + 'gene.png',
             title)
    # family level
    countsVariant = compressCountsFamilyLevel(countsVariant)
    writeCountsToFile(countsVariant,
                      filePrefix + 'family.csv')
    plotDist(countsVariant, sampleName,
             filePrefix + 'family.png',
             title)
