'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import numpy as np
import gc

from math import ceil
from multiprocessing import  Queue
from Bio import SeqIO
from collections import Counter

from abseq.IgRepReporting.igRepPlots import plotDist
from abseq.IgRepAuxiliary.productivityAuxiliary import ProcCounter
from abseq.IgRepAuxiliary.PrimerWorker import PrimerWorker
from abseq.IgRepertoire.igRepUtils import gunzip, calMaxIUPACAlignScores, compressCountsGeneLevel
from abseq.config import MEM_GB


def addPrimerData(cloneAnnot, readFile, format, fr4cut, trim5end,
                  trim3end, actualQstart, end5, end3, end5offset, threads):
    print("Primer specificity analysis has begun ...")
    queryIds = cloneAnnot.index
    seqsPerFile = 100
    _addPrimerColumns(cloneAnnot, end5, end3)
    workers = []
    records = SeqIO.index(gunzip(readFile), format)
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
    except Exception as e:
        print("Something went wrong during the primer specificity analysis!")
        raise e
    finally:
        for w in workers:
            w.terminate()
        records.close()


def _addPrimerColumns(cloneAnnot, end5, end3):
    if end5:
        cloneAnnot['5end'] = np.nan
        cloneAnnot['5endPrimer'] = np.nan
        cloneAnnot['5endIndel'] = np.nan
    if end3:
        cloneAnnot['3end'] = np.nan
        cloneAnnot['3endPrimer'] = np.nan
        cloneAnnot['3endIndel'] = np.nan


def parsePrimerFile(primerFile):
    if primerFile:
        primerSequences = []
        maxPrimerlength = float('-inf')
        primerids = []
        for rec in SeqIO.parse(primerFile, "fasta"):
            maxPrimerlength = max(maxPrimerlength, len(rec.seq))
            primerids.append(rec.id)
            primerSequences.append(str(rec.seq).upper())
        maxScores = calMaxIUPACAlignScores(primerSequences)
        return maxPrimerlength, zip(primerids, primerSequences, maxScores)
    return None, None


def write5EndPrimerStats(name, cloneAnnot, fileprefix, category="All"):
    valid5End = Counter(cloneAnnot['5end'].tolist())
    plotDist(valid5End, name, fileprefix + 'integrity_dist.png',
             title='Integrity of 5`-end Primer Sequence (%s)' % (category),
             proportion=True, rotateLabels=False)
    invalid5Clones = cloneAnnot.index[cloneAnnot['5end'] == 'Indelled'].tolist()
    print("Example of Indelled 5`-end:", invalid5Clones[1:10])
    print("Example of valid 5`-end:", cloneAnnot.index[cloneAnnot['5end'] != 'Indelled'].tolist()[1:10])

    stopcodonInFrameDist = Counter(cloneAnnot['stopcodon'].tolist())
    plotDist(stopcodonInFrameDist, name,
             fileprefix + 'stopcodon_dist.png',
             title='Stop Codons in 5`-End (%s)' % (category),
             proportion=False, rotateLabels=False)

    c1 = Counter(cloneAnnot[cloneAnnot['5end'] == 'Indelled']['5endPrimer'].tolist())
    plotDist(c1, name, fileprefix +
             'indelled_dist.png',
             title='Abundance of Indelled 5`-end Primers (%s)' % (category),
             proportion=False, rotateLabels=False, vertical=False, top=50)
    c = Counter(cloneAnnot[cloneAnnot['5end'] == 'Indelled']['5endIndel'].tolist())
    plotDist(c, name, fileprefix +
             'indel_pos_dist.png',
             title='Abundance of Indel Positions in 5`-end Primers (%s)' % (category),
             proportion=False, rotateLabels=False, vertical=True,
             sortValues=False, top=50)
    primers = set(cloneAnnot['5endPrimer'].tolist())
    # print(c1, primers)
    for primer in primers:
        # print(primer)
        df = cloneAnnot[cloneAnnot['5end'] == 'Indelled']
        df = df[df['5endPrimer'] == primer]
        # print(df.shape)
        germLineDist = compressCountsGeneLevel(Counter(df['vgene'].tolist()))
        plotDist(germLineDist, name, fileprefix + primer +
                 '_igv_dist.png',
                 title='IGV Abundance (%s)' % (category),
                 proportion=False, vertical=False, top=20, rotateLabels=False)
    gc.collect()


def writePrimerStats(end, name, cloneAnnot, fileprefix, category="All"):
    validEnd = Counter(cloneAnnot['{}end'.format(end)].tolist())

    plotDist(validEnd, name, fileprefix + 'integrity_dist.png',
             title='Integrity of {}\'-end Primer Sequence (%s)'.format(end) % (category),
             proportion=True, rotateLabels=False)

    invalidClones = cloneAnnot.index[cloneAnnot['{}end'.format(end)] == 'Indelled'].tolist()

    print("Example of Indelled {}'-end:".format(end), invalidClones[1:10])
    print("Example of valid {}'-end:".format(end),
          cloneAnnot.index[cloneAnnot['{}end'.format(end)] != 'Indelled'].tolist()[1:10])

    stopcodonInFrameDist = Counter(cloneAnnot['stopcodon'].tolist())
    plotDist(stopcodonInFrameDist, name,
             fileprefix + 'stopcodon_dist.png',
             title='Stop Codons in {}\'-End ({})'.format(end, category),
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
    gc.collect()
