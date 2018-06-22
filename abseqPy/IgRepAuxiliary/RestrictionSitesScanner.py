'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import sys

from multiprocessing import Process
from numpy import isnan
from Bio.Seq import Seq

import abseqPy.IgRepAuxiliary.restrictionAuxiliary
from abseqPy.logger import printto, LEVEL


class RestrictionSitesScanner(Process):
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
            except Exception as e:
                printto(self.stream, "An error occurred while processing " + self.name, LEVEL.ERR)
                # raise
                # sys.exit()
                self.resultsQueue.put(None)
                continue
        return
    
    def runSimple(self, nextTask):
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
            record = record[qstart:end]  
            seq = record
            seqRC = str(Seq(seq).reverse_complement())            
            cut = False
            for site in stats["siteHitsCount"].keys():
                hits = abseqPy.IgRepAuxiliary.restrictionAuxiliary.findHits(seq, self.sites[site])
                if len(hits) == 0:                    
                    hits = abseqPy.IgRepAuxiliary.restrictionAuxiliary.findHits(seqRC, self.sites[site])
                if len(hits) > 0:
                    stats["siteHitsCount"][site] += len(hits) 
                    stats["siteHitSeqsCount"][site] += 1                     
                    stats["siteHitsSeqsIDs"][site].append(id_)
                    cut = True                 
            if cut:
                stats["seqsCutByAny"] += 1
        self.procCounter.increment(len(nextTask))                         
        self.resultsQueue.put(stats)
        
    def runDetailed(self, nextTask):
        pass
    
  
def sliceRecord(rec, qsRec):
    qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
    if isnan(qsRec['fr4.end']):
        end = len(rec.seq)
    else:
        end = int(qsRec['fr4.end'])
    return rec[qstart:end]  
