from collections import Counter
from math import ceil
from os.path import exists
import os
from multiprocessing import Queue
from igBlastWorker import analyzeSmallFile, IgBlastWorker
import sys
from Bio import SeqIO
from pandas.core.frame import DataFrame
import gc
from config import MEM_GB
from igRepUtils import splitFastaFile, writeListToFile, writeCountsToFile,\
    compressCountsGeneLevel, compressCountsFamilyLevel
from igRepPlots import plotDist, plotStatsHeatmap


def analyzeIGSeqRead(igRep, fastaFile, bitScore, alignLen, subjStart,
                               seqType='dna'):        
        noWorkers = igRep.threads
        
        if (fastaFile == None):
            return Counter()
        # Estimate the IGV diversity in a library from igblast output 
        print('The IGV abundance of ' + fastaFile.split('/')[-1] + ' is being analyzed ...')
        with open(fastaFile) as f:
            noSeqs = sum(1 for line in f if line.startswith(">"))    
        totalFiles = int(ceil(noSeqs / (igRep.seqsPerFile)))   
#         print(fastaFile, totalFiles, igRep.seqsPerFile)
#         sys.exit()
        if (totalFiles == 1):
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
            (cdrInfo, fileteredIDs) = analyzeSmallFile(newFastFile, igRep.chain, igRep.db,
                                                 bitScore, alignLen, subjStart,
                                                  seqType, noWorkers)
            sys.stdout.flush()
        else:
            # split FASTA file into smaller files 
            ext = '.' + fastaFile.split('/')[-1].split('.')[-1]
            filesDir = igRep.outputDir + "tmp"
            prefix = fastaFile.split('/')[-1].split('.')[0]
            prefix = prefix[prefix.find("_R")+1:prefix.find("_R")+3] + "_" if (prefix.find("_R") != -1) else ""
            splitFastaFile(fastaFile, totalFiles, igRep.seqsPerFile, 
                           filesDir, igRep.primer, prefix, ext)               

            # # Prepare the multiprocessing queues     
            tasks = Queue()    
            outcomes = Queue()   
            exitQueue = Queue()
            if (noWorkers > totalFiles):
                noWorkers = totalFiles          
            cdrInfo = DataFrame()
            fileteredIDs = []
            try:    
                # Initialize workers 
                workers = []        
                for i in range(noWorkers):
                    w = IgBlastWorker(igRep.chain, igRep.db, bitScore, alignLen, subjStart,
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
                    (cdrInfoi, fileteredIDsi) = outcome
                    cdrInfo = cdrInfo.append(cdrInfoi)
                    fileteredIDs += fileteredIDsi
                    sys.stdout.flush()
                    gc.collect()
                    
                    
            except Exception:
                print("Something went wrong during the analysis process!")
                raise
            finally:
                for w in workers:
                    w.terminate()
        #     print("here3")
            # Clean folders to save space
            if (noSeqs > igRep.seqsPerFile and os.path.exists(filesDir + "/part1" + ext)): 
                os.system("rm " + filesDir + "/*" + ext)        
        return (cdrInfo, fileteredIDs)
    

def writeIGVAbundanceToFiles(stats, sampleName, outDir):
        ighvDist = Counter(stats["vgene"].tolist())
        if (len(ighvDist) == 0):
            print("WARNING: No IGV hits were detected.")
            return        
        
        # Write the counts of all IGVs into a text file
        writeCountsToFile(ighvDist, outDir + sampleName + 
                          '_igv_dist_variant_level.csv')
        
        # Group IGVs based on the subfamilies (gene level) and then write into a text file
        ighvDistSub = compressCountsGeneLevel(ighvDist)
#         for k in ighvDist.keys():
#             ksub = k.split('*')[0]
#             ighvDistSub[ksub] = ighvDistSub.get(ksub, 0) + ighvDist[k]
        writeCountsToFile(ighvDistSub, outDir + sampleName + 
                          '_igv_dist_gene_level.csv')
        plotDist(ighvDistSub, sampleName, outDir + sampleName + 
                 '_igv_dist_gene_level.png', rotateLabels=False, vertical=False)
        
        # Group IGVs based on the families and then write into a text file
        ighvDistfam = compressCountsFamilyLevel(ighvDistSub)
#         for k in ighvDistSub.keys():
#             kfam = k.split('-')[0].split('/')[0]
#             ighvDistfam[kfam] = ighvDistfam.get(kfam, 0) + ighvDistSub[k]
        writeCountsToFile(ighvDistfam, outDir + sampleName + 
                           '_igv_dist_family_level.csv')
        
        # Plot the family level distribution
        plotDist(ighvDistfam, sampleName, outDir + sampleName + 
                 '_igv_dist_family_level.png')
        
        # plot alignment length vs %identity
        plotStatsHeatmap(stats, sampleName, ['alignlen', 'identity'],
                         ['Alignment Length', '%Identity'] , outDir + sampleName + 
                  '_align_quality_identity_hm.png')
        # plot alignment length vs bitScore
        plotStatsHeatmap(stats, sampleName, ['alignlen', 'bitscore'],
                         ['Alignment Length', 'bitScore'] , outDir + sampleName + 
                  '_align_quality_bitscore_hm.png')
        # plot query start vs. subject start
        plotStatsHeatmap(stats, sampleName, ['qstart', 'sstart'],
                         ['Query Start', 'Subject Start'] , outDir + sampleName + 
                  '_align_quality_start_hm.png')
        plotStatsHeatmap(stats, sampleName, ['alignlen', 'mismatches'],
                         ['Alignment Length', 'Mismatches'] , outDir + sampleName + 
                  '_align_quality_mismatches_hm.png')
        c = Counter(stats['mismatches'].tolist())
        plotDist(c, sampleName, outDir + sampleName + 
                 '_mismatches_dist.png', title='Number of Mismatches in V gene',
                 proportion=True, rotateLabels=False, top=20) 
        plotStatsHeatmap(stats, sampleName, ['alignlen', 'gaps'],
                         ['Alignment Length', 'Gaps'] , outDir + sampleName + 
                  '_align_quality_gaps_hm.png')
        c = Counter(stats['gaps'].tolist())
        plotDist(c, sampleName, outDir + sampleName + 
                 '_gaps_dist.png', title='Number of Gaps in V gene',
                 proportion=True, rotateLabels=False, top=20) 
    #     print(np.percentile(stats, [0, 100], 0))
    #     summarizeStats(stats, outputDir+sampleName+'_stats_summary.txt')
    
    
