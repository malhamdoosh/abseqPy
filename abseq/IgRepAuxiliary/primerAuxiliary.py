'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import numpy as np
import gc
import sys

from math import ceil
from multiprocessing import Queue
from Bio import SeqIO
from collections import Counter
from pandas import DataFrame

from abseq.IgRepAuxiliary.productivityAuxiliary import ProcCounter
from abseq.IgRepReporting.igRepPlots import plotDist, plotVenn
from abseq.IgRepAuxiliary.PrimerWorker import PrimerWorker
from abseq.IgRepertoire.igRepUtils import gunzip, compressCountsGeneLevel
from abseq.config import MEM_GB



def addPrimerData(cloneAnnot, readFile, format, fr4cut, trim5end,
                  trim3end, actualQstart, end5, end3, end5offset, threads):
    print("Primer specificity analysis has begun ...")
    queryIds = cloneAnnot.index
    seqsPerFile = 100
    _addPrimerColumns(cloneAnnot, end5, end3)
    workers = []
    records = SeqIO.index(gunzip(readFile), format)
    newColumns = ['queryid'] + list(cloneAnnot.columns)
    try:
        print("\t " + format + " index created and refinement started ...")
        noSeqs = len(queryIds)
        totalTasks = int(ceil(noSeqs * 1.0 / seqsPerFile))
        tasks = Queue()
        exitQueue = Queue()
        resultsQueue = Queue()
        procCounter = ProcCounter(noSeqs)
        threads = min(threads, totalTasks)
        if MEM_GB < 16:
            threads = 2
        for _ in range(threads):
            w = PrimerWorker(procCounter, fr4cut, trim5end, trim3end, actualQstart, end5,
                             end3, end5offset, tasks, exitQueue, resultsQueue)
            workers.append(w)
            w.start()
        for i in range(totalTasks):
            ids = queryIds[i * seqsPerFile:(i+1)*seqsPerFile]
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
        print("All workers have completed their tasks successfully.")
        print("Results are being collated from all workers ...")
        cloneAnnotList = _collectPrimerResults(newColumns, resultsQueue, totalTasks, noSeqs)
        print("Results were collated successfully.")
    except Exception as e:
        print("Something went wrong during the primer specificity analysis!")
        raise e
    finally:
        for w in workers:
            w.terminate()
        records.close()

    primerAnnot = DataFrame(cloneAnnotList, columns=newColumns)
    primerAnnot.index = primerAnnot.queryid
    del primerAnnot['queryid']
    return primerAnnot


def _collectPrimerResults(columns, queue, totalTasks, noSeqs):
    processed = 0
    cloneAnnot = []
    while totalTasks:
        result = queue.get()
        totalTasks -= 1
        if result is None:
            continue
        for entry in result:
            # put them as a list (in the ordering specified by 'columns')
            cloneAnnot.append([entry[col] for col in columns])
        processed = len(cloneAnnot)
        if processed % 50000 == 0:
            print("\t{:,}/{:,} records have been collected ... ".format(processed, noSeqs))
            sys.stdout.flush()

    print("\t{:,}/{:,} records have been collected ... ".format(processed, noSeqs))
    return cloneAnnot


def _addPrimerColumns(cloneAnnot, end5, end3):
    if end5:
        cloneAnnot['5end'] = str(np.nan)
        cloneAnnot['5endPrimer'] = str(np.nan)
        cloneAnnot['5endIndel'] = np.nan
    if end3:
        cloneAnnot['3end'] = str(np.nan)
        cloneAnnot['3endPrimer'] = str(np.nan)
        cloneAnnot['3endIndel'] = np.nan


def writePrimerStats(end, name, cloneAnnot, fileprefix, category="All"):

    validEnd = Counter(cloneAnnot['{}end'.format(end)].tolist())

    plotDist(validEnd, name, fileprefix + 'integrity_dist.png',
             title='Integrity of {}\'-end Primer Sequence (%s)'.format(end) % (category),
             proportion=True, rotateLabels=False)

    invalidClones = cloneAnnot.index[cloneAnnot['{}end'.format(end)] == 'Indelled'].tolist()

    print("Example of Indelled {}'-end: {}".format(end, str(invalidClones[1:10])))
    print("Example of valid {}'-end: {}".format(end,
          str(cloneAnnot.index[cloneAnnot['{}end'.format(end)] != 'Indelled'].tolist()[1:10])))

    stopcodonInFrameDist = Counter(cloneAnnot['stopcodon'].tolist())
    plotDist(stopcodonInFrameDist, name,
             fileprefix + 'stopcodon_dist.png',
             title='Stop Codons in sequences by {}\'-End ({})'.format(end, category),
             proportion=False, rotateLabels=False)

    c1 = Counter(cloneAnnot[cloneAnnot['{}end'.format(end)] == 'Indelled']['{}endPrimer'.format(end)].tolist())
    plotDist(c1, name, fileprefix +
             'indelled_dist.png',
             title='Abundance of Indelled {}\'-end Primers ({})'.format(end, category),
             proportion=False, rotateLabels=False, vertical=False, top=50)
    c = Counter(cloneAnnot[cloneAnnot['{}end'.format(end)] == 'Indelled']['{}endIndel'.format(end)].tolist())
    plotDist(c, name, fileprefix +
             'indel_pos_dist.png',
             title='Abundance of Indel Positions in {}\'-end Primers ({})'.format(end, category),
             proportion=False, rotateLabels=False, vertical=True,
             sortValues=False, top=50)
    primers = set(cloneAnnot['{}endPrimer'.format(end)].tolist())
    # print(c1, primers)
    for primer in primers:
        # print(primer)
        df = cloneAnnot[cloneAnnot['{}end'.format(end)] == 'Indelled']
        df = df[df['{}endPrimer'.format(end)] == primer]
        # print(df.shape)
        germLineDist = compressCountsGeneLevel(Counter(df['vgene'].tolist()))
        plotDist(germLineDist, name, fileprefix + primer +
                 '_igv_dist.png',
                 title='IGV Abundance (%s)' % (category),
                 proportion=False, vertical=False, top=20, rotateLabels=False)


def generatePrimerPlots(cloneAnnot, outDir, name, end5, end3):
    nanString = 'NaN'
    # similar with productivity analysis etc ..
    cloneAnnot.fillna(nanString, inplace=True)

    outOfFrameClones = cloneAnnot[cloneAnnot['v-jframe'] == 'Out-of-frame']
    productiveClones = cloneAnnot[(cloneAnnot['v-jframe'] == 'In-frame') & (cloneAnnot['stopcodon'] == 'No')]

    if end5:
        print("5-end analysis of all clones ... ")
        writePrimerStats('5', name, cloneAnnot, outDir + name + '_all_5end_')
        allInvalid5Clones = cloneAnnot.index[cloneAnnot['5end'] == 'Indelled'].tolist()

        print('5-end analysis of out-of-frame clones ... ')
        writePrimerStats('5', name, outOfFrameClones, outDir + name + '_outframe_5end_', 'Out-of-frame')
        outFrameInvalid5Clones = outOfFrameClones.index[outOfFrameClones['5end'] == 'Indelled'].tolist()

        print("5-end analysis of productive clones ... ")
        writePrimerStats('5', name, productiveClones, outDir + name + '_productive_5end_', 'Productive')
        productiveInvalid5Clones = productiveClones.index[productiveClones['5end'] == 'Indelled'].tolist()

    if end3:
        print("3-end analysis of all clones ... ")
        writePrimerStats('3', name, cloneAnnot, outDir + name + '_all_3end_')

        print("3-end analysis of out-of-frame clones ... ")
        writePrimerStats('3', name, outOfFrameClones, outDir + name + '_outframe_3end_', 'Out-of-frame')

        print('3-end analysis of productive clones ... ')
        writePrimerStats('3', name, productiveClones, outDir + name + "_productive_3end_", 'Productive')

        if end5:
            invalid3Clones = cloneAnnot.index[cloneAnnot['3end'] == 'Indelled'].tolist()
            plotVenn({"5'-end": set(allInvalid5Clones), "3'-end": set(invalid3Clones)},
                     outDir + name + '_all_invalid_primers.png')
            del invalid3Clones, allInvalid5Clones

            outFrameInvalid3Clones = outOfFrameClones.index[outOfFrameClones['3end'] == 'Indelled'].tolist()
            plotVenn({"5'-end": set(outFrameInvalid5Clones), "3'-end": set(outFrameInvalid3Clones)},
                     outDir + name + '_outframe_invalid_primers.png')
            del outFrameInvalid3Clones, outFrameInvalid5Clones

            productiveInvalid3Clones = productiveClones.index[productiveClones['3end'] == 'Indelled'].tolist()
            plotVenn({"5'-end": set(productiveInvalid5Clones), "3'-end": set(productiveInvalid3Clones)},
                     outDir + name + "_productive_invalid_primers.png")

    # similar with abundance analysis etc ..
    cloneAnnot.replace(nanString, np.nan, inplace=True)
    gc.collect()
