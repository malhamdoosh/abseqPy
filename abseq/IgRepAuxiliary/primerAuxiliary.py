'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import numpy as np

from math import ceil
from multiprocessing import  Queue
from Bio import SeqIO

from abseq.IgRepAuxiliary.productivityAuxiliary import ProcCounter
from abseq.IgRepAuxiliary.PrimerWorker import PrimerWorker
from abseq.IgRepertoire.igRepUtils import gunzip, calMaxIUPACAlignScores
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


