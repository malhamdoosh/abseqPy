'''
Created on 05/08/2016

@author: monther
'''
from collections import Counter
from math import ceil
import os
from multiprocessing import Queue
from IgRepAuxiliary.IgBlastWorker import analyzeSmallFile, IgBlastWorker
import sys
from Bio import SeqIO
from pandas.core.frame import DataFrame
import gc
from IgRepertoire.igRepUtils import splitFastaFile



        
def annotateIGSeqRead(igRep, fastaFile, seqType='dna'):        
        noWorkers = igRep.threads
        seqsPerFile = igRep.seqsPerFile
        if (fastaFile == None):
            return Counter()
        # Estimate the IGV diversity in a library from igblast output 
        print('The IGV clones of ' + fastaFile.split('/')[-1] + ' are being annotated ...')
        with open(fastaFile) as f:
            noSeqs = sum(1 for line in f if line.startswith(">"))
        totalFiles = int(ceil(noSeqs / seqsPerFile))  
        if totalFiles <  noWorkers:
            seqsPerFile = int(noSeqs * 1.0 / noWorkers) 
            totalFiles = int(ceil(noSeqs * 1.0 / seqsPerFile))
        print("\t{0:,} sequences were found to be distributed into {1:,} files".format(noSeqs,
                                                                             totalFiles))  
#         print(noSeqs, seqsPerFile, totalFiles)
#         sys.exit()
        if igRep.primer > 0:
            recordsAll = SeqIO.to_dict(SeqIO.parse(fastaFile, 'fasta'))
            records = []
            for id in recordsAll:
                rec = recordsAll[id]
                rec.description = ''
                rec.seq = rec.seq[:igRep.primer]
                records.append(rec)
            filesDir = igRep.outputDir + "tmp"
            SeqIO.write(records, filesDir + "/seqs.fasta", 'fasta')
            newFastFile = filesDir + "/seqs.fasta"
        else:
            newFastFile = fastaFile
        if (noWorkers == 1):            
            (cloneAnnot, fileteredIDs) = analyzeSmallFile(newFastFile, igRep.chain, igRep.db,                                                 
                                                  seqType, noWorkers)
            sys.stdout.flush()
        else:
            # split FASTA file into smaller files 
            ext = '.' + fastaFile.split('/')[-1].split('.')[-1]
            filesDir = igRep.outputDir + "tmp"
            prefix = fastaFile.split('/')[-1].split('.')[0]
            prefix = prefix[prefix.find("_R")+1:prefix.find("_R")+3] + "_" if (prefix.find("_R") != -1) else ""
            splitFastaFile(fastaFile, totalFiles, seqsPerFile, 
                           filesDir, prefix, ext)               

            # # Prepare the multiprocessing queues     
            tasks = Queue()    
            outcomes = Queue()   
            exitQueue = Queue()              
            cloneAnnot = DataFrame()
            fileteredIDs = []
            try:    
                # Initialize workers 
                workers = []        
                for i in range(noWorkers):
                    w = IgBlastWorker(igRep.chain, igRep.db,
                                      seqType, int(ceil(noWorkers * 1.0/ totalFiles)))
                    w.tasksQueue = tasks
                    w.resultsQueue = outcomes
                    w.exitQueue = exitQueue      
                    workers.append(w)
                    w.start()       
                    sys.stdout.flush()
                # initialize tasks queue with file names     
                if (noSeqs > igRep.seqsPerFile): 
                    for i in range(totalFiles):
                        tasks.put(filesDir + "/" + prefix + "part"  + `int(i + 1)` + ext)
                else:
                    tasks.put(fastaFile)
                # Add a poison pill for each worker
                for i in range(noWorkers + 10):
                    tasks.put(None)                  
               
                # Wait all process workers to terminate    
                i = 0 
                while i < noWorkers:    
                    m = exitQueue.get()
                    if m == "exit":
                        i += 1
                
                # Collect results
                print("Results are being collated from all workers ...")  
                sys.stdout.flush()
                while totalFiles:
                    outcome = outcomes.get()
                    totalFiles -= 1                    
                    if (outcome is None):
                        continue                    
                    (cloneAnnoti, fileteredIDsi) = outcome
                    cloneAnnot = cloneAnnot.append(cloneAnnoti)
                    fileteredIDs += fileteredIDsi
                    sys.stdout.flush()
                    gc.collect()
                print("\tResults were collated successfully.")
                    
                    
            except Exception:
                print("Something went wrong during the annotation process!")
                raise
            finally:
                for w in workers:
                    w.terminate()
        #     print("here3")
            # Clean folders to save space
            #TODO: remove .fasta and .out files 
            if (noSeqs > igRep.seqsPerFile and 
                os.path.exists(filesDir + "/" + prefix + "part1" + ext)): 
                os.system("rm " + filesDir + "/*" + ext)        
        return (cloneAnnot, fileteredIDs)
    
    
    
    