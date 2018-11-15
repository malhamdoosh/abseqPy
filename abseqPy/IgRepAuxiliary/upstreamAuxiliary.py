'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
import os
import itertools
import gc

from numpy import Inf, random
from Bio import SeqIO
from collections import defaultdict

from abseqPy.IgRepReporting.igRepPlots import plotSeqLenDist, plotSeqLenDistClasses, plotDist
from abseqPy.IgRepertoire.igRepUtils import gunzip, compressCountsFamilyLevel, \
    compressCountsGeneLevel, safeOpen, compressSeqGeneLevel, compressSeqFamilyLevel
from abseqPy.logger import LEVEL, printto
from abseqPy.utilities import requires

_UPSTREAM_SEQ_FILE_SEP = '|'
_VALID_SEQ_FASTA_TEMPLATE = "{}_{}_{:.0f}_{:.0f}_valid_seqs.fasta"
_FAULTY_SEQ_FASTA_TEMPLATE = "{}_{}_{:.0f}_{:.0f}_faulty_trans.fasta"
_STARTCOD_SEQ_FASTA_TEMPLATE = "{}_{}_{:.0f}_{:.0f}_no_atg.fasta"


def plotUpstreamLenDist(upstreamFile, expectLength, name, outDir, stream=None):
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

    :param outDir: string
                output directory

    :param stream:
                logging stream

    :return:
                None
    """
    if not upstreamFile.endswith(".fasta"):
        raise ValueError("Expected {} to be a FASTA file".format(os.path.basename(upstreamFile.rstrip(os.path.sep))))

    fileExt = 'fasta'

    outputFile = os.path.join(outDir, os.path.basename(upstreamFile).replace('.fasta', '_dist.csv'))
    plotSeqLenDist(upstreamFile, name, outputFile, fileExt, stream=stream)

    outputFile = os.path.join(outDir, os.path.basename(upstreamFile).replace('.fasta', '_dist_class.csv'))
    plotSeqLenDistClasses(upstreamFile, name, outputFile, fileExt, stream=stream)

    # if user provided values to upstream (and it's not Inf)
    if expectLength != Inf:
        # plot seqs that do not meet the expectedLength (shorter than expected length)
        outputFile = os.path.join(outDir, os.path.basename(upstreamFile).replace('.fasta', '_dist_short.csv'))
        plotSeqLenDist(upstreamFile, name, outputFile, fileExt, expectLength - 1, stream=stream)

        outputFile = os.path.join(outDir, os.path.basename(upstreamFile).replace('.fasta', '_dist_short_class.csv'))
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
                    record.id = record.id + _UPSTREAM_SEQ_FILE_SEP + qsRec.vgene
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


def collectUpstreamSeqs(upstreamFile, sampleName, expectLength, outResDir, outHdfDir,
                        startCodon=True, type='secsig', plotDist=True, stream=None):
    """
    segregates and plots upstream file sequences. They are segregated as sequences with no start codon,
    faulty sequences (stop codon post translation if type == secsig or X or N nucleotides in the sequence),
    and valid sequences.

    :param upstreamFile: string
                        upstream FASTA file

    :param sampleName: string
                        name of sampel

    :param expectLength: tuple or list
                        index-able of length 2 denoting start and end

    :param outResDir: string
                        name of result output directory

    :param outHdfDir: string
                        name of auxiliary output directory

    :param startCodon: bool
                        whether or not to care about start codons during segregation

    :param type: string
                        either 'secsig' or '5utr'

    :param plotDist: bool
                        whether or not to also save a txt and png file denoting the distribution of segregated sequences

    :param stream: stream
                        debugging stream

    :return: tuple
                        (ighvValidSignals : dict, faultySeqs : dict and noStartCodonSeqs: dict)
    """
    if type not in ['secsig', '5utr']:
        raise ValueError("Unknown parameter type={}, expected one of 'secsig', '5utr'".format(type))

    printto(stream, "\tSequences between {} and {} are being extracted ... ".format(expectLength[0], expectLength[1]))

    START_CODON = "ATG"

    # valid sequences
    ighvSignals = defaultdict(list)
    ighvSignalsCounts = defaultdict(int)

    # no start codons
    ighvSignalsNoATG = defaultdict(list)
    noStartCodonCounts = defaultdict(int)

    # faulty translations
    faultyTrans = defaultdict(list)
    faultyTransCounts = defaultdict(int)

    ignoredSeqs = 0

    records = SeqIO.index(gunzip(upstreamFile), 'fasta')
    for id_ in records:
        rec = records[id_]
        ighv = rec.id.split(_UPSTREAM_SEQ_FILE_SEP)[1]
        seq = rec.seq
        if expectLength[0] <= len(rec) <= expectLength[1]:
            if not startCodon or START_CODON in seq:

                if type == 'secsig':
                    seq = seq[: len(seq) - (len(seq) % 3)].translate(to_stop=False)[1:]

                if 'X' in seq or '*' in seq:
                    faultyTrans[ighv].append(rec)
                    faultyTransCounts[ighv] += 1
                elif 'N' not in rec.seq:
                    ighvSignals[ighv].append(rec)
                    ighvSignalsCounts[ighv] += 1
                else:
                    printto(stream, "Ignored: " + str(rec.seq) + ' ' + str(seq))
                    if type == 'secsig':
                        faultyTrans[ighv].append(rec)
                        faultyTransCounts[ighv] += 1
            elif startCodon:
                # START_CODON not in seq
                ighvSignalsNoATG[ighv].append(rec)
                noStartCodonCounts[ighv] += 1
        else:
            ignoredSeqs += 1

    if ignoredSeqs:
        printto(stream, "\tThere are {} sequences that were ignored because the length of the provided upstream"
                        "sequences were not {} <= length(upstream_seqs) <= {}"
                .format(ignoredSeqs, *expectLength),
                LEVEL.WARN)

    if sum(ighvSignalsCounts.values()):
        flattenRecs = list(itertools.chain.from_iterable(ighvSignals.values()))
        assert len(flattenRecs) == sum(ighvSignalsCounts.values())
        title = 'Valid Secretion Signals' if type == 'secsig' else "Valid 5'-UTRs"
        printto(stream, "\tThere are {} {} within expected "
                        "length ({} to {}) and startCodon={}"
                .format(sum(ighvSignalsCounts.values()), title, expectLength[0], expectLength[1], startCodon),
                LEVEL.INFO)
        validSeqFile = os.path.join(outHdfDir, _VALID_SEQ_FASTA_TEMPLATE
                                    .format(sampleName, type, *expectLength))
        SeqIO.write(flattenRecs, validSeqFile, 'fasta')
        if plotDist:
            writeCountsCategoriesToFile(ighvSignalsCounts,
                                        sampleName,
                                        os.path.join(outResDir, "{}_{}_{:.0f}_{:.0f}_valid_"
                                                     .format(sampleName, type, expectLength[0], expectLength[1])),
                                        title)
    if sum(faultyTransCounts.values()):
        flattenRecs = (itertools.chain.from_iterable(faultyTrans.values()))
        assert len(flattenRecs) == sum(faultyTransCounts.values())
        faultySeqFile = os.path.join(outHdfDir, _FAULTY_SEQ_FASTA_TEMPLATE
                                     .format(sampleName, type, *expectLength))
        SeqIO.write(flattenRecs, faultySeqFile, 'fasta')
        if plotDist:
            writeCountsCategoriesToFile(faultyTransCounts,
                                        sampleName,
                                        os.path.join(outResDir, "{}_{}_{:.0f}_{:.0f}_faulty_"
                                                     .format(sampleName, type, *expectLength)),
                                        'Faulty Translations')
        printto(stream, "\tTotal faulty secretion signals is {} (excluded)".format(len(flattenRecs)), LEVEL.INFO)
        for i in random.choice(range(len(flattenRecs)), min(5, len(flattenRecs)), replace=False):
            sequence = flattenRecs[i].seq
            printto(stream,
                    "\t{}\n\tTranslated:{}"
                    .format(sequence, sequence[: len(sequence) - (len(sequence) % 3)].translate()))

    if sum(noStartCodonCounts.values()):
        flattenRecs = list(itertools.chain.from_iterable(ighvSignalsNoATG.values()))
        assert len(flattenRecs) == sum(noStartCodonCounts.values())
        noStartCodonFile = os.path.join(outHdfDir, _STARTCOD_SEQ_FASTA_TEMPLATE
                                        .format(sampleName, type, *expectLength))
        SeqIO.write(flattenRecs, noStartCodonFile, 'fasta')
        if plotDist:
            writeCountsCategoriesToFile(noStartCodonCounts,
                                        sampleName,
                                        os.path.join(outResDir, "{}_{}_{:.0f}_{:.0f}_no_atg_"
                                                     .format(sampleName, type, *expectLength)),
                                        "Upstream sequences without start codon")
        printto(stream,
                "\tThere is no ATG codon in {} sequences (excluded)".format(len(flattenRecs)),
                LEVEL.INFO)
        for i in random.choice(range(len(flattenRecs)), min(5, len(flattenRecs)), replace=False):
            printto(stream, "\t{}".format(flattenRecs[i].seq))

    # the output of each ighv key's value should be a list of strings, not SeqRecord object
    for k in ighvSignals:
        ighvSignals[k] = map(lambda x: str(x.seq), ighvSignals[k])
    for k in faultyTrans:
        faultyTrans[k] = map(lambda x: str(x.seq), faultyTrans[k])
    for k in ighvSignalsNoATG:
        ighvSignalsNoATG[k] = map(lambda x: str(x.seq), ighvSignalsNoATG[k])

    return ighvSignals, faultyTrans, ighvSignalsNoATG


# generateMotifs will be skipped silently without an exception if TAMO is not found
@requires("TAMO", fatal=False)
def findUpstreamMotifs(upstreamFile, sampleName, outHdfDir, outResDir, expectLength, level,
                       startCodon=True, type='secsig', clusterMotifs=False, threads=2, stream=None):
    """
    finds and visualizes motifs from the sequences provided in upstreamFile

    :param upstreamFile: string
                    path to FASTA file containing upstream sequences

    :param sampleName: string
                    name to refer the sample as

    :param outHdfDir: string
                    path to aux directory

    :param outResDir: string
                    path to result directory

    :param expectLength: tuple or list
                    index-able of length 2 denoting start and end.
                    If start == end, this implies that the analysis should
                    be conducted ONLY on sequences with length == start == end, the rest are ignored.

    :param level: string
                    one of 'gene', 'family' or 'variant'

    :param startCodon: bool
                    whether or not to segregate sequences with start codon

    :param type: string
                    one of upstream analysis types: '5utr' or 'secsig'

    :param clusterMotifs: bool
                    whether or not to cluster sequences using TAMO

    :param threads: int
                    number of threads to use

    :param stream: stream
                    logging stream
    :return: None
    """
    from abseqPy.IgRepAuxiliary.seqUtils import generateMotifs

    if level == 'variant':
        # single argument identity function
        compressor = lambda signals: signals
    elif level == 'gene':
        compressor = compressSeqGeneLevel
    elif level == 'family':
        compressor = compressSeqFamilyLevel
    else:
        raise ValueError("Unknown level {} requested, accepted values are family, gene, or variant".format(level))

    if type not in ['secsig', '5utr']:
        raise ValueError("Unknown parameter type={}, expected one of 'secsig', '5utr'".format(type))

    # output files always have this format: <sampleName>_<type>_<exp[0]>_<exp[1]>_*
    OUTPUT_FILE_PACKET = (sampleName, type, expectLength[0], expectLength[1])

    # only analyze motifs of secretion signals that have exactly length == expectLength[0] == expectLength[1]
    EXACT_LENGTH = expectLength[0] == expectLength[1]

    validSeqFile = os.path.join(outHdfDir, _VALID_SEQ_FASTA_TEMPLATE.format(*OUTPUT_FILE_PACKET))
    faultySeqFile = os.path.join(outHdfDir, _FAULTY_SEQ_FASTA_TEMPLATE.format(*OUTPUT_FILE_PACKET))
    noStartCodonFile = os.path.join(outHdfDir, _STARTCOD_SEQ_FASTA_TEMPLATE.format(*OUTPUT_FILE_PACKET))

    allFiles = [validSeqFile, faultySeqFile, noStartCodonFile]

    if all(map(lambda x: os.path.exists(x), allFiles)):
        printto(stream,
                "Sequences were already analyzed at {}, loading from files instead ... " + ' '.join(allFiles),
                LEVEL.WARN)

        ighvSignals, faultySeq, noStartCodonSeq = _loadIGVSeqsFromFasta(validSeqFile),\
                                                  _loadIGVSeqsFromFasta(faultySeqFile),\
                                                  _loadIGVSeqsFromFasta(noStartCodonFile)
    else:
        printto(stream, "Sequences are being analyzed ... ")
        ighvSignals, faultySeq, noStartCodonSeq = collectUpstreamSeqs(upstreamFile, sampleName, expectLength,
                                                                      outResDir, outHdfDir, startCodon, type,
                                                                      stream=stream)

    ighvSignals = compressor(ighvSignals)
    generateMotifs(ighvSignals,
                   align=(expectLength[0] < expectLength[1]),
                   outputPrefix=os.path.join(outResDir,
                                             ("{}_{}_{:.0f}_{:.0f}_dna_" + level).format(*OUTPUT_FILE_PACKET)),
                   clusterMotifs=clusterMotifs,
                   threads=threads,
                   stream=stream)

    if EXACT_LENGTH and type == 'secsig':
        faultySeq = compressor(faultySeq)
        generateMotifs(faultySeq,
                       align=True,
                       outputPrefix=os.path.join(outResDir,
                                                 ("{}_{}_{:.0f}_{:.0f}_faulty_" + level).format(*OUTPUT_FILE_PACKET)),
                       transSeq=False,
                       extendAlphabet=True,
                       clusterMotifs=clusterMotifs,
                       threads=threads,
                       stream=stream)
        noStartCodonSeq = compressor(noStartCodonSeq)
        generateMotifs(noStartCodonSeq,
                       align=True,
                       outputPrefix=os.path.join(outResDir,
                                                 ("{}_{}_{:.0f}_{:.0f}_untranslated_" + level).format(
                                                     *OUTPUT_FILE_PACKET)),
                       transSeq=False,
                       extendAlphabet=True,
                       clusterMotifs=clusterMotifs,
                       threads=threads,
                       stream=stream)
        generateMotifs(ighvSignals,
                       align=False,
                       outputPrefix=os.path.join(outResDir,
                                                 ("{}_{}_{:.0f}_{:.0f}_protein_" + level).format(*OUTPUT_FILE_PACKET)),
                       transSeq=True,
                       clusterMotifs=clusterMotifs,
                       threads=threads,
                       stream=stream)


def writeCountsCategoriesToFile(countsVariant, sampleName, filePrefix, title=''):
    # gene level
    countsVariant = compressCountsGeneLevel(countsVariant)
    plotDist(countsVariant, sampleName, filePrefix + 'gene.csv', title)
    # family level
    countsVariant = compressCountsFamilyLevel(countsVariant)
    plotDist(countsVariant, sampleName, filePrefix + 'family.csv', title)


def _loadIGVSeqsFromFasta(filename):
    ighvSeqs = defaultdict(list)
    with safeOpen(filename) as fp:
        for rec in SeqIO.parse(fp, 'fasta'):
            ighv = rec.id.split('|')[1].strip()
            ighvSeqs[ighv].append(str(rec.seq))
    return ighvSeqs
