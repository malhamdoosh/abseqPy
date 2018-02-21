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

import abseq.IgRepAuxiliary.restrictionAuxiliary


class RestrictionSitesScanner(Process):
    '''
    Restriction Sites Scanner
    '''
    def __init__(self, records, cloneAnnot, procCounter, sites, simpleScan = True):
        '''
        Constructor
        '''
        super(RestrictionSitesScanner, self).__init__()
        self.records = records
        self.cloneAnnot = cloneAnnot
        self.procCounter = procCounter
        self.sites = sites
        self.simpleScan = simpleScan
        self.tasksQueue = None
        self.exitQueue = None
        self.resultsQueue = None   
        
    def run(self):
        print(self.name + " process is now ready to start a new job ...")
        sys.stdout.flush() 
        while True:            
            nextTask = self.tasksQueue.get()
            #t = time.time()
            if (nextTask is None):
                print(self.name + " process has stopped." )
                self.exitQueue.put("exit")
                break
            #print("Waiting for a job took: {0:f}".format(time.time() - t))    
            try:
                #t = time.time()
                if self.simpleScan:
                    self.runSimple(nextTask)
                else:
                    self.runDetailed(nextTask)
                #print("Running a job took: {0:f}".format(time.time() - t))    
            except Exception as e:
                print("An error occurred while processing " + self.name)
                print(e)
                #raise
                sys.exit()
                self.resultsQueue.put(None)
                continue
        return
    
    def runSimple(self, nextTask):
        stats = abseq.IgRepAuxiliary.restrictionAuxiliary.initSimpleRSAStats(self.sites)
        stats['total'] = len(nextTask)      
        for id in nextTask:
            record = self.records[id]
            qsRec = self.cloneAnnot.loc[id].to_dict()
            qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
            if (isnan(qsRec['fr4.end'])):
                end = len(record)
            else:
                end = int(qsRec['fr4.end'])
            record = record[qstart:end]  
            seq = record
            seqRC = str(Seq(seq).reverse_complement())            
            cut = False
            for site in stats["siteHitsCount"].keys():
                hits = abseq.IgRepAuxiliary.restrictionAuxiliary.findHits(seq, self.sites[site])
                if len(hits) == 0:                    
                    hits = abseq.IgRepAuxiliary.restrictionAuxiliary.findHits(seqRC, self.sites[site])
                if len(hits) > 0:
                    stats["siteHitsCount"][site] += len(hits) 
                    stats["siteHitSeqsCount"][site] += 1                     
                    stats["siteHitsSeqsIDs"][site].append(id)   
                    cut = True                 
            if cut:
                stats["seqsCutByAny"] += 1
        self.procCounter.increment(len(nextTask))                         
        self.resultsQueue.put(stats)
        
    def runDetailed(self, nextTask):
        pass
    
  
def sliceRecord((rec, qsRec)):
    qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
    if (isnan(qsRec['fr4.end'])):
        end = len(rec.seq)
    else:
        end = int(qsRec['fr4.end'])
    return rec[qstart:end]  
