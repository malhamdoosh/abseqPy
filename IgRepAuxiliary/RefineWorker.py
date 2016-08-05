'''
Created on 05/08/2016

@author: monther
'''
from multiprocessing import Process


class RefineWorker(Process):
    def __init__(self, chain, igBlastDB, 
                 seqType, threads):
        super(RefineWorker, self).__init__() 
        self.chain = chain     
        self.threads = threads
        self.tasksQueue = None
        self.resultsQueue = None
        self.exitQueue = None
    
    def run(self):
        while True:            
            nextTask = self.tasksQueue.get()
            # poison pill check            
            if (nextTask is None):
                print("process has stopped ... " + self.name)
                self.exitQueue.put("exit")
                break
            try:
                result = None
#                 analyzeSmallFile(nextTask, self.chain, self.igBlastDB,                                                                                                      
#                                                      self.seqType, self.threads)                        
#                 print("process has completed analysis... " + self.name) 
                self.resultsQueue.put(result)            
            except Exception as e:
                print("An error occurred while processing " + nextTask.split('/')[-1])
                print(e)
                self.resultsQueue.put(None)
#                 raise
#                 sys.exit()
                continue                       
#             print("process has completed a run... " + self.name) 
        return