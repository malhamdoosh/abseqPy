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
import numpy as np

def annotateIGSeqRead(igRep, fastaFile, seqType='dna'):        
        noWorkers = igRep.threads
        
        if (fastaFile == None):
            return Counter()
        # Estimate the IGV diversity in a library from igblast output 
        print('The IGV clones of ' + fastaFile.split('/')[-1] + ' are being annotated ...')
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
            #TODO: remove .fasta and .out files 
            if (noSeqs > igRep.seqsPerFile and 
                os.path.exists(filesDir + "/" + prefix + "part1" + ext)): 
                os.system("rm " + filesDir + "/*" + ext)        
        return (cdrInfo, fileteredIDs)
    

def writeVAbundanceToFiles(stats, sampleName, outDir):
    igvDist = Counter(stats["vgene"].tolist())
    if (len(igvDist) == 0):
        print("WARNING: No IGV hits were detected.")
        return        
    
    # Write the counts of all IGVs into a text file
    writeCountsToFile(igvDist, outDir + sampleName + 
                      '_igv_dist_variant_level.csv')
    
    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igvDistSub = compressCountsGeneLevel(igvDist)
#         for k in igvDist.keys():
#             ksub = k.split('*')[0]
#             igvDistSub[ksub] = igvDistSub.get(ksub, 0) + igvDist[k]
    writeCountsToFile(igvDistSub, outDir + sampleName + 
                      '_igv_dist_gene_level.csv')
    plotDist(igvDistSub, sampleName, outDir + sampleName + 
             '_igv_dist_gene_level.png', rotateLabels=False, vertical=False)
    
    # Group IGVs based on the families and then write into a text file
    igvDistfam = compressCountsFamilyLevel(igvDistSub)
#         for k in igvDistSub.keys():
#             kfam = k.split('-')[0].split('/')[0]
#             igvDistfam[kfam] = igvDistfam.get(kfam, 0) + igvDistSub[k]
    writeCountsToFile(igvDistfam, outDir + sampleName + 
                       '_igv_dist_family_level.csv')
    
    # Plot the family level distribution
    plotDist(igvDistfam, sampleName, outDir + sampleName + 
             '_igv_dist_family_level.png')
    
    # plot alignment length vs %identity
    plotStatsHeatmap(stats, sampleName, ['alignlen', 'identity'],
                     ['Alignment Length', '%Identity'] , outDir + sampleName + 
              '_igv_align_quality_identity_hm.png')
    # plot alignment length vs bitScore
    plotStatsHeatmap(stats, sampleName, ['alignlen', 'bitscore'],
                     ['Alignment Length', 'bitScore'] , outDir + sampleName + 
              '_igv_align_quality_bitscore_hm.png')
    # plot query start vs. subject start
    plotStatsHeatmap(stats, sampleName, ['vqstart', 'vstart'],
                     ['Query Start', 'Subject Start'] , outDir + sampleName + 
              '_igv_align_quality_start_hm.png')
    plotStatsHeatmap(stats, sampleName, ['alignlen', 'vmismatches'],
                     ['Alignment Length', 'Mismatches'] , outDir + sampleName + 
              '_igv_align_quality_mismatches_hm.png')
    c = Counter(stats['vmismatches'].tolist())
    plotDist(c, sampleName, outDir + sampleName + 
             '_igv_mismatches_dist.png', title='Number of Mismatches in V gene',
             proportion=True, rotateLabels=False, top=20) 
    plotStatsHeatmap(stats, sampleName, ['alignlen', 'vgaps'],
                     ['Alignment Length', 'Gaps'] , outDir + sampleName + 
              '_igv_align_quality_gaps_hm.png')
    c = Counter(stats['vgaps'].tolist())
    plotDist(c, sampleName, outDir + sampleName + 
             '_igv_gaps_dist.png', title='Number of Gaps in V gene',
             proportion=True, rotateLabels=False, top=20) 
    #     print(np.percentile(stats, [0, 100], 0))
    #     summarizeStats(stats, outputDir+sampleName+'_stats_summary.txt')


def writeJAbundanceToFiles(stats, sampleName, outDir):
    igjDist = Counter(stats["jgene"].tolist())
    igjDist = {str(k) : igjDist[k] for k in igjDist}
    if (len(igjDist) == 0):
        print("WARNING: No IGJ hits were detected.")
        return        
    
    # Write the counts of all IGVs into a text file
    writeCountsToFile(igjDist, outDir + sampleName + 
                      '_igj_dist_variant_level.csv')
    
    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igjDistSub = compressCountsGeneLevel(igjDist)
#     writeCountsToFile(igjDistSub, outDir + sampleName + 
#                       '_igj_dist_gene_level.csv')
#     plotDist(igjDistSub, sampleName, outDir + sampleName + 
#              '_igj_dist_gene_level.png', rotateLabels=False, vertical=False)
#     
    # Group IGVs based on the families and then write into a text file
    igjDistfam = compressCountsFamilyLevel(igjDistSub)
    writeCountsToFile(igjDistfam, outDir + sampleName + 
                       '_igj_dist_family_level.csv')    
    # Plot the family level distribution
    plotDist(igjDistfam, sampleName, outDir + sampleName + 
             '_igj_dist_family_level.png',
             title = 'IGJ Abundance in Sample ' + sampleName )


def writeDAbundanceToFiles(stats, sampleName, outDir):
    igdDist = Counter(stats["dgene"].tolist())
    igdDist = Counter({str(k) : igdDist[k] for k in igdDist})
    if (len(igdDist) == 0):
        print("WARNING: No IGD hits were detected.")
        return        
    
    # Write the counts of all IGVs into a text file
    writeCountsToFile(igdDist, outDir + sampleName + 
                      '_igd_dist_variant_level.csv')
    
    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igdDistSub = compressCountsGeneLevel(igdDist)
    writeCountsToFile(igdDistSub, outDir + sampleName + 
                      '_igd_dist_gene_level.csv')
    plotDist(igdDistSub, sampleName, outDir + sampleName + 
             '_igd_dist_gene_level.png', rotateLabels=False, vertical=False,
             title = 'IGD Abundance in Sample ' + sampleName )
    
    # Group IGVs based on the families and then write into a text file
    igdDistfam = compressCountsFamilyLevel(igdDistSub)
    writeCountsToFile(igdDistfam, outDir + sampleName + 
                       '_igd_dist_family_level.csv')    
    # Plot the family level distribution
    plotDist(igdDistfam, sampleName, outDir + sampleName + 
             '_igd_dist_family_level.png',
             title = 'IGD Abundance in Sample ' + sampleName)


def writeAbundanceToFiles(stats, sampleName, outDir, chain = "hv"):
        writeVAbundanceToFiles(stats, sampleName, outDir)
        writeJAbundanceToFiles(stats, sampleName, outDir)
        if (chain == "hv"):
            writeDAbundanceToFiles(stats, sampleName, outDir)
    
    
