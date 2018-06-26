'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''


from multiprocessing import Process
from numpy import isnan
from Bio.Seq import Seq

import abseqPy.IgRepAuxiliary.restrictionAuxiliary
from abseqPy.logger import printto, LEVEL


class RestrictionSitesScanner(Process):
    # todo: remove cloneAnnot. using (xthreads) expensive memory multiplier
    def __init__(self, records, cloneAnnot, procCounter, sites, simpleScan=True, stream=None):
        super(RestrictionSitesScanner, self).__init__()
        self.records = records
        self.cloneAnnot = cloneAnnot
        self.procCounter = procCounter
        self.sites = sites
        self.simpleScan = simpleScan
        self.tasksQueue = None
        self.exitQueue = None
        self.resultsQueue = None
        self.stream = stream

    def run(self):
        printto(self.stream, self.name + " process is now ready to start a new job ...")
        while True:
            nextTask = self.tasksQueue.get()
            if nextTask is None:
                printto(self.stream, self.name + " process has stopped.")
                self.exitQueue.put("exit")
                break
            try:
                if self.simpleScan:
                    self.runSimple(nextTask)
                else:
                    self.runDetailed(nextTask)
            except Exception:
                printto(self.stream, "An error occurred while processing " + self.name, LEVEL.ERR)
                self.resultsQueue.put(None)
                continue
        return
    
    def runSimple(self, nextTask):
        """
        Runs Restriction sites simple analysis

        :param nextTask: iterable of sequence ids that should exist in self.records
        :return: None
        """
        stats = abseqPy.IgRepAuxiliary.restrictionAuxiliary.initSimpleRSAStats(self.sites)
        stats['total'] = len(nextTask)
        for id_ in nextTask:
            # record = raw sequence (taken from m.dict())
            record = self.records[id_]
            qsRec = self.cloneAnnot.loc[id_].to_dict()
            qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
            if isnan(qsRec['fr4.end']):
                end = len(record)
            else:
                end = int(qsRec['fr4.end'])
            seq = record[qstart:end]
            seqRC = str(Seq(seq).reverse_complement())
            cut = False
            for site, siteRegex in self.sites.items():
                hits = abseqPy.IgRepAuxiliary.restrictionAuxiliary.findHits(seq, siteRegex)
                if len(hits) == 0:                    
                    hits = abseqPy.IgRepAuxiliary.restrictionAuxiliary.findHits(seqRC, siteRegex)
                if len(hits) > 0:
                    # how many times has this site found a match on this sequence (seq / seqRC)
                    stats["siteHitsCount"][site] += len(hits)
                    # how many times has this site found a match on *a* sequence
                    # (yes, this is a "duplicate" field of siteHitsSeqsIDs, we could've taken the length of
                    # siteHitsSeqsIDs, it would be equal to this)
                    stats["siteHitSeqsCount"][site] += 1
                    # add the ids of sequences where this site has a match for
                    stats["siteHitsSeqsIDs"][site].append(id_)
                    cut = True                 
            if cut:
                # total number of sequences that are cut by *at least* one site
                stats["seqsCutByAny"] += 1
        self.procCounter.increment(len(nextTask))                         
        self.resultsQueue.put(stats)
        
    def runDetailed(self, nextTask):
        # todo: move code here
        pass
    
  
def sliceRecord(rec, qsRec):
    qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
    if isnan(qsRec['fr4.end']):
        end = len(rec.seq)
    else:
        end = int(qsRec['fr4.end'])
    return rec[qstart:end]  
