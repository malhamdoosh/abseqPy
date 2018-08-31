'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import os
import gc
import sys
import numpy as np

from math import ceil
from multiprocessing import Queue
from Bio import SeqIO
from collections import Counter
from pandas import DataFrame

from abseqPy.IgRepAuxiliary.productivityAuxiliary import ProcCounter
from abseqPy.IgRepReporting.igRepPlots import plotDist, plotVenn
from abseqPy.IgRepAuxiliary.PrimerWorker import PrimerWorker
from abseqPy.IgRepertoire.igRepUtils import gunzip, compressCountsGeneLevel
from abseqPy.utilities import hasLargeMem
from abseqPy.logger import printto, LEVEL


def addPrimerData(cloneAnnot, readFile, format, fr4cut, trim5end,
                  trim3end, actualQstart, end5, end3, end5offset, threads, stream=None):
    printto(stream, "Primer specificity analysis has begun ...")
    queryIds = cloneAnnot.index
    seqsPerFile = 100
    _addPrimerColumns(cloneAnnot, end5, end3)
    workers = []
    records = SeqIO.index(gunzip(readFile), format)
    newColumns = ['queryid'] + list(cloneAnnot.columns)
    try:
        printto(stream, "\t " + format + " index created and primer analysis started ...")
        noSeqs = len(queryIds)
        totalTasks = int(ceil(noSeqs * 1.0 / seqsPerFile))
        tasks = Queue()
        exitQueue = Queue()
        resultsQueue = Queue()
        procCounter = ProcCounter(noSeqs, stream=stream)
        threads = min(threads, totalTasks)
        if not hasLargeMem():
            threads = 2
        for _ in range(threads):
            w = PrimerWorker(procCounter, fr4cut, trim5end, trim3end, actualQstart, end5,
                             end3, end5offset, tasks, exitQueue, resultsQueue, stream=stream)
            workers.append(w)
            w.start()
        for i in range(totalTasks):
            ids = queryIds[i * seqsPerFile:(i + 1) * seqsPerFile]
            recs = map(lambda x: records[x], ids)
            qsRecs = map(lambda x: cloneAnnot.loc[x].to_dict(), ids)
            tasks.put((recs, qsRecs))

        # poison pills
        for _ in range(threads + 10):
            tasks.put(None)

        i = 0
        while i < threads:
            m = exitQueue.get()
            i += (m == 'exit')
        printto(stream, "All workers have completed their tasks successfully.")
        printto(stream, "Results are being collated from all workers ...")
        cloneAnnotList = _collectPrimerResults(newColumns, resultsQueue, totalTasks, noSeqs, stream=stream)
        printto(stream, "Results were collated successfully.")
    except Exception as e:
        printto(stream, "Something went wrong during the primer specificity analysis!", LEVEL.EXCEPT)
        raise e
    finally:
        for w in workers:
            w.terminate()
        records.close()

    primerAnnot = DataFrame(cloneAnnotList, columns=newColumns)
    primerAnnot.set_index('queryid', drop=True, inplace=True)
    return primerAnnot


def _collectPrimerResults(columns, queue, totalTasks, noSeqs, stream=None):
    processed = 0
    cloneAnnot = []
    totalUnexpected5 = totalUnexpected3 = 0
    while totalTasks:
        result = queue.get()
        totalTasks -= 1
        if result is None:
            continue
        for entry, unexpected5, unexpected3 in result:
            totalUnexpected5 += unexpected5
            totalUnexpected3 += unexpected3
            # put them as a list (in the ordering specified by 'columns')
            cloneAnnot.append([entry[col] for col in columns])
        processed = len(cloneAnnot)
        if processed % 50000 == 0:
            printto(stream, "\t{:,}/{:,} records have been collected ... ".format(processed, noSeqs))
            sys.stdout.flush()

    printto(stream, "\t{:,}/{:,} records have been collected ... ".format(processed, noSeqs))
    printto(stream, "\tThere were {} unexpected 5' alignment and {} unexpected 3' alignment"
            .format(totalUnexpected5, totalUnexpected3), LEVEL.WARN)
    return cloneAnnot


def _addPrimerColumns(cloneAnnot, end5, end3):
    if end5:
        cloneAnnot['5endPrimer'] = str(np.nan)
        cloneAnnot['5endMismatchIndex'] = np.nan
        cloneAnnot['5endIndelIndex'] = np.nan
    if end3:
        cloneAnnot['3endPrimer'] = str(np.nan)
        cloneAnnot['3endMismatchIndex'] = np.nan
        cloneAnnot['3endIndelIndex'] = np.nan


def writePrimerStats(end, name, cloneAnnot, fileprefix, category="All", stream=None):
    NA = str(np.nan)
    PRIMER = str(end) + 'endPrimer'
    MISMATCH = str(end) + 'endMismatchIndex'
    INDEL = str(end) + 'endIndelIndex'

    known = cloneAnnot[cloneAnnot[PRIMER] != NA]
    integrity = {
        'Unknown': (len(cloneAnnot) - len(known)),
        'Indelled': sum(known[INDEL] != 0),
        'Mismatched': sum(known[MISMATCH] != 0),
        'Intact': len(known[(known[INDEL] == 0) & (known[MISMATCH] == 0)])
    }

    plotDist(integrity, name, fileprefix + 'integrity_dist.csv',
             title='Integrity of {}\'-end Primer Sequence (%s)'.format(end) % (category),
             proportion=True, rotateLabels=False)

    invalidClones = known.index[known[INDEL] != 0].tolist()
    valid = known.index[known[INDEL] == 0].tolist()
    printto(stream, "Example of Indelled {}'-end: {}".format(end, str(invalidClones[1:10])), LEVEL.INFO)
    printto(stream, "Example of non-indelled {}'-end: {}".format(end, str(valid[1:10])), LEVEL.INFO)

    c1 = Counter(known[known[INDEL] != 0][PRIMER].tolist())
    plotDist(c1, name, fileprefix +
             'indelled_dist.csv',
             title='Abundance of Indelled {}\'-end Primers ({})'.format(end, category),
             proportion=False, rotateLabels=False, vertical=False, top=50)

    c = Counter(known[known[INDEL] != 0][INDEL].tolist())
    plotDist(c, name, fileprefix +
             'indel_pos_dist.csv',
             title='Abundance of Indel Positions in {}\'-end Primers ({})'.format(end, category),
             proportion=False, rotateLabels=False, vertical=True,
             sortValues=False, top=50)

    primers = set(known[PRIMER].tolist())

    for primer in primers:
        # get only ighv abundance of indelled primers
        df = known[known[INDEL] != 0]
        df = df[df[PRIMER] == primer]

        germLineDist = compressCountsGeneLevel(Counter(df['vgene'].tolist()))
        plotDist(germLineDist, name, fileprefix + primer +
                 '_igv_dist.csv',
                 title='IGV Abundance of indelled {} ({})'.format(primer, category),
                 proportion=False, vertical=False, top=20, rotateLabels=False)


def generatePrimerPlots(cloneAnnot, outDir, name, end5, end3, stream=None):
    nanString = 'NaN'
    # similar with productivity analysis etc ..
    cloneAnnot.fillna(nanString, inplace=True)
    NA = str(np.nan)
    PRIMER5 = '5endPrimer'
    PRIMER3 = '3endPrimer'
    INDEL5 = '5endIndelIndex'
    INDEL3 = '3endIndelIndex'

    outOfFrameClones = cloneAnnot[cloneAnnot['v-jframe'] == 'Out-of-frame']
    productiveClones = cloneAnnot[(cloneAnnot['v-jframe'] == 'In-frame') & (cloneAnnot['stopcodon'] == 'No')]

    if end5:
        printto(stream, "5-end analysis of all clones ... ")
        writePrimerStats('5', name, cloneAnnot, os.path.join(outDir, name + '_all_5end_'), stream=stream)
        allInvalid5Clones = cloneAnnot.index[
            ((cloneAnnot[PRIMER5] != NA) & (cloneAnnot[INDEL5] != 0))
        ].tolist()

        printto(stream, '5-end analysis of out-of-frame clones ... ')
        writePrimerStats('5', name, outOfFrameClones, os.path.join(outDir, name + '_outframe_5end_'), 'Out-of-frame',
                         stream=stream)
        outFrameInvalid5Clones = outOfFrameClones.index[
            ((outOfFrameClones[PRIMER5] != NA) & (outOfFrameClones[INDEL5] != 0))
        ].tolist()

        printto(stream, "5-end analysis of productive clones ... ")
        writePrimerStats('5', name, productiveClones, os.path.join(outDir, name + '_productive_5end_'), 'Productive',
                         stream=stream)
        productiveInvalid5Clones = productiveClones.index[
            ((productiveClones[PRIMER5] != NA) & (productiveClones[INDEL5] != 0))
        ].tolist()

    if end3:
        printto(stream, "3-end analysis of all clones ... ")
        writePrimerStats('3', name, cloneAnnot, os.path.join(outDir, name + '_all_3end_'), stream=stream)

        printto(stream, "3-end analysis of out-of-frame clones ... ")
        writePrimerStats('3', name, outOfFrameClones, os.path.join(outDir, name + '_outframe_3end_'), 'Out-of-frame',
                         stream=stream)

        printto(stream, '3-end analysis of productive clones ... ')
        writePrimerStats('3', name, productiveClones, os.path.join(outDir, name + "_productive_3end_"), 'Productive',
                         stream=stream)

        if end5:
            invalid3Clones = cloneAnnot.index[
                ((cloneAnnot[PRIMER3] != NA) & (cloneAnnot[INDEL3] != 0))
            ].tolist()
            plotVenn({"5'-end": set(allInvalid5Clones), "3'-end": set(invalid3Clones)},
                     os.path.join(outDir, name + '_all_invalid_primers.png'),
                     "Intersection of indelled 5' and 3' sequences (All)", stream=stream)
            del invalid3Clones, allInvalid5Clones

            outFrameInvalid3Clones = outOfFrameClones.index[
                ((outOfFrameClones[PRIMER3] != NA) & (outOfFrameClones[INDEL3] != 0))
            ].tolist()
            plotVenn({"5'-end": set(outFrameInvalid5Clones), "3'-end": set(outFrameInvalid3Clones)},
                     os.path.join(outDir, name + '_outframe_invalid_primers.png'),
                     "Intersection of indelled 5' and 3' sequences (Out-of-frame)", stream=stream)
            del outFrameInvalid3Clones, outFrameInvalid5Clones

            productiveInvalid3Clones = productiveClones.index[
                ((productiveClones[PRIMER3] != NA) & (productiveClones[INDEL3] != 0))
            ].tolist()
            plotVenn({"5'-end": set(productiveInvalid5Clones), "3'-end": set(productiveInvalid3Clones)},
                     os.path.join(outDir, name + "_productive_invalid_primers.png"),
                     "Intersection of indelled 5' and 3' sequences (productive)", stream=stream)
    # similar with abundance analysis etc ..
    cloneAnnot.replace(nanString, np.nan, inplace=True)
    gc.collect()
