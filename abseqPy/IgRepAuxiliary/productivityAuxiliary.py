'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from __future__ import division
import gc
import sys
import os

from Bio import SeqIO
from collections import defaultdict, Counter
from pandas.core.frame import DataFrame
from numpy import random, isnan
from multiprocessing import Queue, Value, Lock
from math import ceil

from abseqPy.IgRepAuxiliary.RefineWorker import RefineWorker
from abseqPy.IgRepertoire.igRepUtils import gunzip
from abseqPy.IgRepAuxiliary.IgBlastWorker import getAnnotationFields
from abseqPy.logger import LEVEL, printto
from abseqPy.utilities import hasLargeMem


def loadRefineFlagInfo():
    refineFlagMsgs = {
        'fr1NotAtBegin': "{:,} clones have FR1 start not equal to query start (Excluded)",

        'endsWithStopCodon': "{:,} clones contain a stop codon ",
        'updatedStopCodon': "The stopcodon flag was updated for {:,} clones ",

        'updatedInFrame': "The v-j frame rearrangement status has been corrected for {:,} clones ",
        'updatedInFrameNA': "{:,} clones have undefined in-frame status",
        'updatedInFrameConc': "{:,} clones show discordance between the query and v gene starts",
        'updatedInFrameNo3or4': "{:,} clones have no CDR3 or FR4",
        'updatedInFrame3x': "{:,} clones are not multiple of 3 ",
        'updatedInFrameIndel': "{:,} clones have indels in one of the FRs or CDRs",

        'CDR3dna': "The CDR3 of {:,} clones was determined using DNA consensus",

        'partitioning': "{:,} clones were partitioned incorrectly.",

        'FR4PredictedError': "{:,} clones have no FR4 end because the consensus region cannot be identified.",
        # 'FR4AbruptEnd': "{:,} clones' J gene ends abruptly. Sequence ends before end of J gene",        # not in use
        'FR4cutEarly': "{:,} clones have --trim3 sequence(s) match earlier than expected. Matched"
                       " before J germline ends, expected after.",
        'FR4Endless': "{:,} clones do not align with provided trim3 sequences",
        'fr4NotAsExpected': "The FR4 of {:,} clones do not start as expected",
        'noFR4': "{:,} clones do not have FR4 ",
        'filterFRLength': "{:,} clones have different FR1, FR2, FR3, or FR4 length compared to the majority of the "
                          "sequences with similar V/J germline gene (excluded)"
    }
    return list(refineFlagMsgs.keys()), refineFlagMsgs


def refineClonesAnnotation(outDir, sampleName, cloneAnnotOriginal, readFile, format,
                           actualQstart, chain, fr4cut,
                           trim5End, trim3End,
                           seqsPerFile, threads, stream=None):
    printto(stream, "Clone annotation and in-frame prediction are being refined ...")
    seqsPerFile = 100
    cloneAnnot = cloneAnnotOriginal.copy()
    queryIds = cloneAnnot.index
    (refineFlagNames, refineFlagMsgs) = loadRefineFlagInfo()
    records = None
    workers = None
    try:
        # process clones from the FASTA/FASTQ file

        # NOTE:
        # if the readFile is gzipped, we need to unzip it in the same directory before passing into
        # SeqIO.index because it doesn't accept gzipped nor opened files
        records = SeqIO.index(gunzip(readFile), format)
        printto(stream, "\t " + format + " index created and refinement started ...")
        # Parallel implementation of the refinement
        noSeqs = len(queryIds)
        totalTasks = int(ceil(noSeqs / seqsPerFile))
        tasks = Queue()
        exitQueue = Queue()
        resultsQueue = Queue()
        procCounter = ProcCounter(noSeqs, stream=stream)
        threads = min(threads, totalTasks)
        # Initialize workers
        workers = []
        for i in range(threads):
            w = RefineWorker(procCounter, chain, actualQstart, fr4cut,
                             trim5End, trim3End, refineFlagNames, stream=stream)
            w.tasksQueue = tasks
            w.exitQueue = exitQueue
            w.resultsQueue = resultsQueue
            workers.append(w)
            w.start()
            sys.stdout.flush()
            # adding jobs to the tasks queue with subsets of query IDs
        assert (totalTasks >= 1)
        for i in range(totalTasks):
            ids = queryIds[i * seqsPerFile:(i + 1) * seqsPerFile]
            recs = map(lambda x: records[x], ids)
            qsRecs = map(lambda x: cloneAnnot.loc[x].to_dict(), ids)
            tasks.put((recs, qsRecs))
        # Add a poison pill for each worker
        for i in range(threads + 10):
            tasks.put(None)
            # Wait all process workers to terminate                
        i = 0
        while i < threads:
            m = exitQueue.get()
            if m == "exit":
                i += 1
        printto(stream, "All workers have completed their tasks successfully.")
        # Collect results
        printto(stream, "Results are being collated from all workers ...")
        # End of parallel implementation
        sys.stdout.flush()

        # invoking the result collection method
        cloneAnnotList, transSeqs, flags, frameworkLengths = collectRefineResults(resultsQueue, totalTasks,
                                                                                  noSeqs, refineFlagNames,
                                                                                  stream=stream)
        printto(stream, "\tResults were collated successfully.")

        # mark each clone as filtered=yes if its framework len is not the most common among the same V/J germline gene
        printto(stream, "Filtering clones according to framework lengths ... ")
        # display the tallies of all frameworks based on V/J germline gene to log file
        for gene in frameworkLengths:
            printto(stream, "{}:".format(gene), LEVEL.INFO)
            for region, counts in frameworkLengths[gene].items():
                printto(stream, "\t{}: {}".format(region.upper(), str(counts)), LEVEL.INFO)

        # flag them if they're filtered based on FR len
        annotationFields = getAnnotationFields(chain)
        for clone in cloneAnnotList:
            flags['filterFRLength'] += markClones(clone, frameworkLengths, annotationFields)

        filtered = len(flags['filterFRLength'])
        printto(stream, "\t{:.2%} ({}/{}) clones were marked as filtered-out using Framework region 1, 2, 3 "
                        "and 4 lengths".format(filtered / cloneAnnot.shape[0], filtered, cloneAnnot.shape[0]))

        # print refine flags
        printRefineFlags(flags, records, refineFlagNames, refineFlagMsgs, stream=stream)
        printto(stream, "Flagged sequences are being written to an output file ... ")
        writeRefineFlags(flags, records, refineFlagNames, refineFlagMsgs,
                         outDir, sampleName)
    except Exception as e:
        printto(stream, "Something went wrong during the refinement process!", LEVEL.EXCEPT)
        raise e
    finally:
        if workers:
            for w in workers:
                w.terminate()
        if records:
            records.close()

    # Create new data frame of clone annotation
    # add new column for filtering based on FR region
    newColumns = getAnnotationFields(chain) + ['filtered']
    cloneAnnot = DataFrame(cloneAnnotList, columns=newColumns)
    cloneAnnot.set_index('queryid', inplace=True, drop=True)
    gc.collect()

    # Create data frame of FR and CDR sequences
    cols = ['queryid', 'germline', 'fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4']
    cloneSeqs = DataFrame(transSeqs, columns=cols)
    for col in cols:
        cloneSeqs.loc[:, col] = cloneSeqs[col].map(str)
    cloneSeqs.set_index('queryid', inplace=True, drop=True)

    return cloneAnnot, cloneSeqs


def markClones(clone, frameworkLengths, annotationFields):
    id_ = []
    if clone[annotationFields.index('v-jframe')] == 'In-frame':
        if _isExpectedFRLength(annotationFields, frameworkLengths, clone):
            clone.append('No')
        else:
            clone.append('Yes')
            id_ = [clone[annotationFields.index('queryid')]]
    else:
        # not in frame, by default will not be considered by diversity anyway
        clone.append('Yes')
    return id_


def collectRefineResults(resultsQueue, totalTasks, noSeqs, refineFlagNames, stream=None):
    total = 0
    cloneAnnot = []
    transSeqs = []
    frameworkLengths = defaultdict(_defaultCounter)

    flags = {}
    for f in refineFlagNames:
        flags[f] = []

    while totalTasks:
        result = resultsQueue.get()
        totalTasks -= 1
        if result is None:
            continue

        qsRecsOrdered, seqs, flagsi, recordLengths = result

        # convert dict to Counter object
        for keys, regions in recordLengths.items():
            for region in regions:
                frameworkLengths[keys][region] += Counter(recordLengths[keys][region])

        # update relevant annotation fields
        cloneAnnot += qsRecsOrdered
        transSeqs += seqs

        # update flags 
        for f in refineFlagNames:
            flags[f] += flagsi[f]
        total += len(qsRecsOrdered)

        if total % 50000 == 0:
            printto(stream, '\t{}/{} records have been collected ... '.format(total, noSeqs))

    printto(stream, '\t{}/{} records have been collected ... '.format(total, noSeqs))
    return cloneAnnot, transSeqs, flags, frameworkLengths


def printRefineFlags(flags, records, refineFlagNames, refineFlagMsgs, stream=None):
    # print statistics and a few of the flagged clones
    for f in refineFlagNames:
        if len(flags[f]) > 0:
            printto(stream, refineFlagMsgs[f].format(len(flags[f])), LEVEL.INFO)
            examples = random.choice(range(len(flags[f])), min(3, len(flags[f])), replace=False)
            for i in examples:
                printto(stream, ">" + flags[f][i], LEVEL.INFO)
                printto(stream, str(records[flags[f][i]].seq), LEVEL.INFO)


def writeRefineFlags(flags, records, refineFlagNames, refineFlagMsgs, outDir, sampleName):
    # 8gb buffer size if system has large enough memory (-1 implies system buffer size)
    with open(os.path.join(outDir, sampleName + "_refinement_flagged.txt"), 'w',
              buffering=int(1 << 23) if hasLargeMem() else -1) as flaggedFp, \
            open(os.path.join(outDir, sampleName + "_refinement_flagged_summary.txt"), "w") as summaryFp:
        for f in refineFlagNames:
            if len(flags[f]) > 0:
                summaryFp.write(refineFlagMsgs[f].format(len(flags[f])) + "\n")
                flaggedFp.write("# " + refineFlagMsgs[f].format(len(flags[f])) + "\n")
                for i in range(len(flags[f])):
                    flaggedFp.write(">" + flags[f][i] + "\n")
                    flaggedFp.write(str(records[flags[f][i]].seq) + "\n")
                flaggedFp.write("\n")


class ProcCounter(object):
    def __init__(self, noSeqs, initval=0, desc="records", stream=None):
        self.val = Value('i', initval)
        self.desc = desc
        self.noSeqs = noSeqs
        self.lock = Lock()
        self.stream = stream

    def increment(self, val=1):
        with self.lock:
            self.val.value += val
            if self.val.value % 50000 == 0 or self.noSeqs == self.val.value:
                printto(self.stream, '\t{}/{} {} have been processed ... '
                        .format(self.val.value, self.noSeqs, self.desc))

    def value(self):
        with self.lock:
            return self.val.value


def _isExpectedFRLength(annotationFields, germlineFrameworkLength, qsRec):
    vgerm = qsRec[annotationFields.index('vgene')].split('*')[0]
    jgerm = qsRec[annotationFields.index('jgene')].split('*')[0]
    for region in ('fr1', 'fr2', 'fr3', 'fr4'):
        frLength = qsRec[annotationFields.index(region + '.end')] - \
                   qsRec[annotationFields.index(region + '.start')] + 1
        germ = vgerm if region != 'fr4' else jgerm
        if not isnan(frLength) and germlineFrameworkLength[germ][region].most_common(1)[0][0] != frLength:
            return False
    return True


def _defaultCounter():
    """
    to be pick-able, this function cannot be lambda
    :return: equivalent to defaultdict(Counter)
    """
    return defaultdict(Counter)



