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
from config import MEM_GB, FR4_CONSENSUS
from igRepUtils import splitFastaFile, writeListToFile, writeCountsToFile,\
    compressCountsGeneLevel, compressCountsFamilyLevel, extractProteinFrag,\
    findBestAlignment
from igRepPlots import plotDist, plotStatsHeatmap
import numpy as np
from numpy import isnan, random
import traceback

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
            (cloneAnnot, fileteredIDs) = analyzeSmallFile(newFastFile, igRep.chain, igRep.db,                                                 
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
        return (cloneAnnot, fileteredIDs)
    

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
    plotDist(igjDist, sampleName, outDir + sampleName + 
                  '_igj_dist_variant_level.png', rotateLabels=False, vertical=False)
    
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

def refineCloneAnnotation(qsRec, record, actualQstart, fr4cut, chain, flags):
    try:
        seqs = [record.id, qsRec['vgene']]    
        if (qsRec['strand'] == "reversed"):
            record.seq = record.seq.reverse_complement()
        # grab the beginning of the VH clone  
        if actualQstart > -1:
            qstart = actualQstart # zero-based
        else:                                  
            qstart = int(qsRec['vqstart'] - qsRec['vstart'])  # zero-based
        if  qstart < 0:
            qstart = 0                   
        vh = record.seq[qstart:]                       
        # check whether the VH sequence can be translated successfully
        if len(vh) % 3 != 0:
            vh = vh[:-1 * (len(vh) % 3)]#                   
        protein = str(vh.translate())
        # check whether the start of the V gene is the same as the start of FR1
        if qsRec['vqstart'] != qsRec['fr1.start']:
            flags['fr1NotAtBegin'] += [record.id]
        # Extract protein sequence of FR1             
        seqs.append(extractProteinFrag(protein, qstart, qsRec['fr1.end'], qstart))
        # Extract protein sequence of CDR1
        seqs.append(extractProteinFrag(protein, qsRec['cdr1.start'], qsRec['cdr1.end'], qstart))
        # Extract protein sequence of FR2
        seqs.append(extractProteinFrag(protein, qsRec['fr2.start'], qsRec['fr2.end'], qstart))
        # Extract protein sequence of CDR2
        seqs.append(extractProteinFrag(protein, qsRec['cdr2.start'], qsRec['cdr2.end'], qstart))
        # Extract protein sequence of FR3
        seqs.append(extractProteinFrag(protein, qsRec['fr3.start'], qsRec['fr3.end'], qstart))
        # Identification of FR4 so that CDR3 can be defined 
        if isnan(qsRec['fr4.end']):
            fr4start, fr4end = findBestAlignment(
                        extractProteinFrag(protein, qsRec['fr3.end'] + 1,
                     - 1, qstart, trimAtStop=False), FR4_CONSENSUS[chain])#                     
            if (fr4start != -1 and fr4end != -1 and fr4end > fr4start):
                qsRec['fr4.start'] = (fr4start - 1) * 3 + qsRec['fr3.end'] + 1
                ## Check whether to cut the Ig sequence after FR4 or not
                if not fr4cut: 
                    qsRec['fr4.end'] = len(record.seq)  # fr4end * 3 + qsRec['fr3.end']     
                else:
                    qsRec['fr4.end'] = fr4end * 3 + qsRec['fr3.end']                                   
                # CDR3
                qsRec['cdr3.start'] = qsRec['fr3.end'] + 1
                qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
            else:
                #TODO: check this case 
                qsRec['cdr3.start'] = qsRec['fr3.end'] + 1
                qsRec['cdr3.end'] = qsRec['jqend']
        # Extract protein sequence of CDR3 and FR4     
        seqs.append(extractProteinFrag(protein, qsRec['cdr3.start'], qsRec['cdr3.end'], qstart))
        seqs.append(extractProteinFrag(protein, qsRec['fr4.start'], qsRec['fr4.end'], qstart))
        # check whether FR and CDR sequences were extracted correctly
        tmp = ''.join(seqs[2:])
        if (tmp not in protein):
            raise Exception("ERROR at CDR refinement: ", protein, seqs)   
        protein = tmp
                       
        if ('*' in protein):
            flags['endsWithStopCodon'] += [record.id]                     
            ### update the StopCodon value if it was set to No
            if qsRec['stopcodon'] == 'No':
                flags['updatedStopCodon'] += [record.id]
                qsRec['stopcodon'] = 'Yes' 
       
        if (seqs[-1][:len(FR4_CONSENSUS[chain])] != FR4_CONSENSUS[chain]):
            flags['fr4NotAsExpected'] += [record.id]
        if (seqs[-1] is None):
            flags['noFR4'] += [record.id]
        
        #TODO: update the annotation fields with the new calculated values
        qsRec['fr1.start'] = qstart+1                
        gaps = abs(qsRec['vqstart'] - qsRec['vstart']) - qstart                
        mismatches = qsRec['vstart'] - 1
        if (qsRec['vstart'] > qsRec['vqstart'] and gaps > 0):
            mismatches -= gaps
        # Only update gaps if the actual query start position is known 
        if gaps > 0:     
            qsRec['fr1.gaps'] += gaps                    
            qsRec['vgaps'] += gaps
        # if igblast ignores mismatches at the begining ==> update
        if (mismatches > 0):
            qsRec['fr1.mismatches'] += mismatches
            qsRec['vmismatches']  += mismatches
            qsRec['vstart']  -= mismatches
            qsRec['vqstart']  -= mismatches                  
    except Exception as e:                
        print("ERROR: exception in the CDR Annotation Refinement")
        print(protein, record.id, str(vh), qsRec)
        traceback.print_exc(file=sys.stdout)
        raise e
    return seqs

refineFlagNames = ['fr1NotAtBegin', 'endsWithStopCodon', 
             'fr4NotAsExpected', 'updatedStopCodon', 'noFR4',
             'updatedInFrame' , 'updatedInFrameNA', 'updatedInFrameConc',
             'updatedInFrameNo3or4', 'updatedInFrame3x' ,
             'updatedInFrameIndel']
refineFlagMsgs = {}
refineFlagMsgs['fr1NotAtBegin'] = "{:,} clones have FR1 start not equal to query start (Excluded)" 
refineFlagMsgs['endsWithStopCodon'] = "{:,} clones contain a stop codon "
refineFlagMsgs['updatedStopCodon'] = "The stopcodon flag was updated for {:,} clones "
refineFlagMsgs['fr4NotAsExpected'] = "{:,} clones do not have an expected FR4 clones "
refineFlagMsgs['noFR4'] = "{:,} clones do not have FR4 "
refineFlagMsgs['updatedInFrame'] = "The v-j frame rearrangement status has been corrected for {:,} clones "
refineFlagMsgs['updatedInFrameNA'] = "{:,} clones have undefined in-frame status"
refineFlagMsgs['updatedInFrameConc'] = "{:,} clones show discordance between the query and v gene starts"
refineFlagMsgs['updatedInFrameNo3or4'] = "{:,} clones have no CDR3 or FR4"
refineFlagMsgs['updatedInFrame3x'] = "{:,} clones are not multiple of 3 "
refineFlagMsgs['updatedInFrameIndel'] = "{:,} clones have indels in one of the FRs or CDRs"

def refineClonesAnnotation(cloneAnnotOriginal, readFile, format, 
                            actualQstart, chain, fr4cut, seqsPerFile):
        print("Clone annotation and in-frame prediction are being refined ...")
        cloneAnnot = cloneAnnotOriginal.copy()        
        queryIds = cloneAnnot.index
        transSeqs = []
        flags = {}
        for f in refineFlagNames:
            flags[f] = [] 
        procSeqs = 0
        # process clones from the FASTA/FASTQ file
        if (MEM_GB > 20):
            records = SeqIO.to_dict(SeqIO.parse(readFile, format))
        else:
            records = SeqIO.index(readFile, format)
        for id in queryIds:            
            # retrieve the clone record from the CloneAnnot dataframe
            record = records[id] 
            qsRec = cloneAnnot.loc[record.id].to_dict()
            seqs = refineCloneAnnotation(qsRec, record, 
                                          actualQstart, fr4cut, chain, flags)   
            # out-of-frame clones are excluded
            if qsRec['v-jframe'] != 'Out-of-frame':
                refineInFramePrediction(qsRec, record, actualQstart, flags)
            # update relevant annotation fields
            cloneAnnot.set_value(id, 'fr1.gaps', qsRec['fr1.gaps'] )                    
            cloneAnnot.set_value(id, 'vgaps', qsRec['vgaps'])                
            cloneAnnot.set_value(id, 'fr1.mismatches', qsRec['fr1.mismatches'])
            cloneAnnot.set_value(id, 'vmismatches', qsRec['vmismatches'])
            cloneAnnot.set_value(id, 'vstart', qsRec['vstart'] )
            cloneAnnot.set_value(id, 'vqstart', qsRec['vqstart'] )
            cloneAnnot.set_value(id, 'cdr3.start', qsRec['cdr3.start'])
            cloneAnnot.set_value(id, 'cdr3.end', qsRec['cdr3.end'])
            cloneAnnot.set_value(id, 'fr4.start', qsRec['fr4.start'])
            cloneAnnot.set_value(id, 'fr4.end' , qsRec['fr4.end'])  
            cloneAnnot.set_value(id, 'stopcodon', qsRec['stopcodon'])
            cloneAnnot.set_value(id, 'fr1.start', qsRec['fr1.start'])   
            cloneAnnot.set_value(id, 'v-jframe', qsRec['v-jframe'])
            # append the FR and CDR protein clones
            transSeqs.append(seqs)
            procSeqs += 1
            if procSeqs % seqsPerFile == 0:
                print('\t %d/%d clones have been processed ... ' % (procSeqs, len(queryIds)))
                sys.stdout.flush()
        print('%d/%d clones have been processed ... ' % (procSeqs, len(queryIds)))
        # print refine flags 
        printRefineFlags(flags, chain)
        flags = None
        # Create data frame of FR and CDR clones
        cloneSeqs = DataFrame(transSeqs, columns=['queryid', 'germline', 'fr1', 'cdr1', 'fr2', 'cdr2',
                                                  'fr3', 'cdr3', 'fr4'])
        cloneSeqs.index = cloneSeqs.queryid
        del cloneSeqs['queryid']
        return (cloneAnnot, cloneSeqs)
        
def printRefineFlags(flags, chain):
    ## print statistics and a few of the flagged clones 
    for f in refineFlagNames:
        if (len(flags[f]) > 0):
            print(refineFlagMsgs[f].format(len(flags[f])))
            examples = random.choice(range(len(flags[f])), 10)
            for i in examples:
                print(flags[f][i])      
             
def refineInFramePrediction(qsRec, record, actualQstart, flags):
        inframe = True    
        # check the the v-jframe value is not NA       
        if (qsRec['v-jframe'] == 'N/A' or
            (type(qsRec['v-jframe']) != type('str') and isnan(qsRec['v-jframe']))):
            flags['updatedInFrameNA'] += [record.id]
            inframe = False
        # the query clone is not in concordance with the start of the germline gene
        qstart = qsRec['vqstart'] - qsRec['vstart'] + 1 # 1-based 
        if (inframe and (qstart < 1 or 
            (actualQstart != -1 and (qstart - 1 - actualQstart) % 3 != 0))):            
            inframe = False     
            flags['updatedInFrameConc'] += [record.id]           
        # if no CDR3 or FR4 ==> Out-of-frame
        if (inframe and (isnan(qsRec['fr4.start']) or 
            isnan(qsRec['fr4.end'])  or isnan(qsRec['cdr3.start']) or 
           qsRec['cdr3.start'] >= qsRec['cdr3.end'])):
            inframe = False
            flags['updatedInFrameNo3or4'] += [record.id]
        # doesn't start/end properly .. not multiple of 3
        if (inframe and ((qsRec['fr4.end'] - qstart + 1) % 3 != 0 or 
            (actualQstart != -1 and 
             ((qsRec['fr4.end'] - actualQstart) % 3 != 0)))):
            inframe = False
            flags['updatedInFrame3x'] += [record.id]
        # indels (gaps) in FRs or CDRs cause frame-shift ==> out-of-frame
        if (inframe and ( 
            qsRec['fr1.gaps'] % 3 != 0 or 
            qsRec['fr2.gaps'] % 3 != 0 or
            qsRec['fr3.gaps'] % 3 != 0 or 
            qsRec['cdr1.gaps'] % 3 != 0 or 
            qsRec['cdr2.gaps'] % 3 != 0 or
            qsRec['cdr3.gaps'] % 3 != 0)):
            inframe = False            
            flags['updatedInFrameIndel'] += [record.id]
        # FR1 start is not aligned with query start 
#             if (not inframe or qsRec['vqstart'] != qsRec['fr1.start']):
#                 inframe = False
        if not inframe:
            qsRec['v-jframe'] =  'Out-of-frame'
            flags['updatedInFrame'] += [record.id]
        
        
# def refineCloneAnnotation(cloneAnnotOriginal, ):
#         print("CDR3 annotation is being refined ...")
#         cloneAnnot = cloneAnnotOriginal.copy()
#         # loading the 5` and 3` primers and calculate maximum alignment scores
#         if igRep.end5:
#             end5Seqs = [(rec.id, str(rec.seq), len(rec.seq)) for rec in SeqIO.parse(igRep.end5, "fasta")]
#             L5 = max(map(lambda x:x[2], end5Seqs))
#             ids = map(lambda x: x[0], end5Seqs)
#             end5Seqs = map(lambda x: x[1].upper()  , end5Seqs)
#             maxScores = calMaxIUPACAlignScores(end5Seqs)
#             end5Seqs = zip(ids, end5Seqs, maxScores)
#             valid5End = {}
#             primer5End = {}
#             indel5End = {}
#         if igRep.end3:
#             end3Seqs = [(rec.id, str(rec.seq), len(rec.seq)) for rec in SeqIO.parse(igRep.end3, "fasta")]
#             L3 = max(map(lambda x:x[2], end3Seqs))
#             ids = map(lambda x: x[0], end3Seqs)
#             end3Seqs = map(lambda x: x[1].upper()  , end3Seqs)
#             maxScores = calMaxIUPACAlignScores(end3Seqs)
#             end3Seqs = zip(ids, end3Seqs, maxScores)
#             valid3End = {}
#             primer3End = {}
#             indel3End = {}
#         queryIds = cloneAnnot.index
#         transSeqs = []
#         fr1NotAtBegin = []
#         endsWithStopCodon = []
#         fr4NotAsExpected = []
#         updatedStopCodon = 0
#         stopCodonPos = {}
#         noFR4 = []
#         procSeqs = 0
#         protein = ''
#         vh = ''
#         sys.stdout.flush()
#         # process clones from the FASTA file
#         if (MEM_GB > 20):
#             records = SeqIO.to_dict(SeqIO.parse(igRep.readFile1, igRep.format))
#         else:
#             records = SeqIO.index(igRep.readFile1, igRep.format)
#         for id in queryIds:            
#             record = records[id]
#             try:    
#                 # retrive the clone record from the CDRInfo file
#                 qsRec = cloneAnnot.loc[record.id].to_dict()
#                 seqs = [record.id, qsRec['vgene']]                
# #                 if qstart <= igRep.actualQstart:
# #                     continue
#                 if (qsRec['strand'] == "reversed"):
#                     record.seq = record.seq.reverse_complement()
#                 # grab the beginning of the VH clone  
#                 if igRep.actualQstart > -1:
#                     qstart = igRep.actualQstart # zero-based
#                 else:                                  
#                     qstart = int(qsRec['vqstart'] - qsRec['vstart'])  # zero-based
#                 if  qstart < 0:
#                     qstart = 0                   
#                 vh = record.seq[qstart:]                       
#                 # check whether the VH clone can be translated successfully
#                 if len(vh) % 3 != 0:
#                     vh = vh[:-1 * (len(vh) % 3)]#                   
#                 protein = str(vh.translate())
#                 
#                 if qsRec['vqstart'] != qsRec['fr1.start']:
#                     fr1NotAtBegin += [record.id]
# 
#                 # FR1             
#                 seqs.append(extractProteinFrag(protein, qstart,
#                                                qsRec['fr1.end'], qstart))
#                 # CDR1
#                 seqs.append(extractProteinFrag(protein, qsRec['cdr1.start'],
#                                                qsRec['cdr1.end'], qstart))
#                 # FR2
#                 seqs.append(extractProteinFrag(protein, qsRec['fr2.start'],
#                                                qsRec['fr2.end'], qstart))
#                 # CDR2
#                 seqs.append(extractProteinFrag(protein, qsRec['cdr2.start'],
#                                                qsRec['cdr2.end'], qstart))
#                 # FR3
#                 seqs.append(extractProteinFrag(protein, qsRec['fr3.start'],
#                                                qsRec['fr3.end'], qstart))
#                 # Identification of FR4 so that CDR3 can be defined 
#                 if isnan(qsRec['fr4.end']):
#                     fr4start, fr4end = findBestAlignment(
#                                 extractProteinFrag(protein, qsRec['fr3.end'] + 1,
#                              - 1, qstart, trimAtStop=False), FR4_CONSENSUS)#                     
#                     if (fr4start != -1 and fr4end != -1 and fr4end > fr4start):
#                         qsRec['fr4.start'] = (fr4start - 1) * 3 + qsRec['fr3.end'] + 1
#                         if not igRep.fr4cut: 
#                             qsRec['fr4.end'] = len(record.seq)  # fr4end * 3 + qsRec['fr3.end']     
#                         else:
#                             qsRec['fr4.end'] = fr4end * 3 + qsRec['fr3.end']                                   
#                         # CDR3
#                         qsRec['cdr3.start'] = qsRec['fr3.end'] + 1
#                         qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
#                     else:
#                         qsRec['cdr3.start'] = qsRec['fr3.end'] + 1
#                         qsRec['cdr3.end'] = qsRec['jqend']
#                     
#                 seqs.append(extractProteinFrag(protein, qsRec['cdr3.start'],
#                                                    qsRec['cdr3.end'], qstart))
#                 seqs.append(extractProteinFrag(protein, qsRec['fr4.start'],
#                                                qsRec['fr4.end'], qstart))
#                 ## Check whether to cut the Ig clone after FR4 or not
#                 if igRep.fr4cut:
#                     try:                        
#                         protein = ''.join(seqs[2:])
#                     except:
#                         pass
#                     try:
#                         vh = record.seq[qstart:int(qsRec['fr4.end'])]
#                     except:
#                         pass
#                 if ('*' in protein):
# #                     print(protein)
#                     endsWithStopCodon += [record.id] 
#                     # check the location of the stop codon
#                     # (5-end primer, in the middle, 3-end primer)
#                     stopCodonPos[record.id] = []
#                     if '*' in protein[:6] :
#                         stopCodonPos[record.id].append("Yes")
#                     else:
#                         stopCodonPos[record.id].append("No")
#                     if '*' in protein[-6:]:
#                         stopCodonPos[record.id].append("Yes")
#                     else:
#                         stopCodonPos[record.id].append("No")
#                     if '*'in protein[6:-6]:
#                         stopCodonPos[record.id].append("Yes")
#                     else:
#                         stopCodonPos[record.id].append("No")
#                     ### update the StopCodon value if it was set to No
#                     if qsRec['stopcodon'] == 'No':
#                         updatedStopCodon += 1
#                         cloneAnnot.set_value(record.id, 'stopcodon', 'Yes') 
#                 else:
#                     stopCodonPos[record.id] = ["No", "No", "No"]
#                 # check if the primer clones match the 5`-end and 3`-end
#                 
#                 if igRep.end5 and qsRec.get('5end', None) is None:
#                     if igRep.end5offset == 0:
#                         prim = str(vh[:L5])
#                     else:
#                         if qstart + igRep.end5offset >= 0:
#                             prim = str(record.seq[qstart + igRep.end5offset: qstart + L5 + igRep.end5offset])
#                         else:
#                             prim = str(record.seq[: qstart + L5 + igRep.end5offset])
#                     (id, tag, indelPos) = findBestMatchedPattern(prim, end5Seqs)
#                     valid5End[record.id] = tag
#                     primer5End[record.id] = id                    
#                     indel5End[record.id] = indelPos
#                 if igRep.end3 and qsRec.get('3end', None) is None:
#                     (id, tag, indelPos) = findBestMatchedPattern(str(vh[-1*L3:]), end3Seqs) 
#                     valid3End[record.id] = tag
#                     primer3End[record.id] = id                    
#                     indel3End[record.id] = indelPos
#                 if (seqs[-1] != FR4_CONSENSUS):
#                     fr4NotAsExpected += [record.id]
#                 if (seqs[-1] is None):
#                     noFR4 += [record.id]
#                 transSeqs.append(seqs)
#                 #TODO: update the annotation fields with the new calculated values
#                 cloneAnnot.set_value(record.id, 'fr1.start', qstart+1)                
#                 gaps = abs(qsRec['vqstart'] - qsRec['vstart']) - qstart                
#                 mismatches = qsRec['vstart'] - 1
#                 if (qsRec['vstart'] > qsRec['vqstart']):
#                     mismatches -= gaps
#                 # Only update gaps if the actual query start position is known 
#                 if gaps > 0:     
#                     cloneAnnot.set_value(record.id, 'fr1.gaps', qsRec['fr1.gaps'] + gaps)                    
#                     cloneAnnot.set_value(record.id, 'vgaps', qsRec['vgaps'] + gaps)
#                 # if igblast ignores mismatches at the begining ==> update
#                 if (mismatches > 0):
#                     cloneAnnot.set_value(record.id, 'fr1.mismatches', qsRec['fr1.mismatches']  + mismatches)
#                     cloneAnnot.set_value(record.id, 'vmismatches', qsRec['vmismatches']  + mismatches)
#                     cloneAnnot.set_value(record.id, 'vstart', qsRec['vstart']  - mismatches)
#                     cloneAnnot.set_value(record.id, 'vqstart', qsRec['vqstart']  - mismatches)
#                 cloneAnnot.set_value(record.id, 'cdr3.start', qsRec['cdr3.start'])
#                 cloneAnnot.set_value(record.id, 'cdr3.end', qsRec['cdr3.end'])
#                 cloneAnnot.set_value(record.id, 'fr4.start', qsRec['fr4.start'])
#                 cloneAnnot.set_value(record.id, 'fr4.end' , qsRec['fr4.end'])
#                 procSeqs += 1
#                 if procSeqs % igRep.seqsPerFile == 0:
#                     print('%d/%d clones have been processed ... ' % (procSeqs, len(queryIds)))
#                     sys.stdout.flush()
#             except Exception as e:                
#                 print("ERROR: exception in the CDR Annotation Refinement")
#                 print(protein, record.id, str(vh), qsRec)
#                 traceback.print_exc(file=sys.stdout)
#                 raise e
# 
#         print('%d/%d clones have been processed ... ' % (procSeqs, len(queryIds)))
#         # Expand the CDRInfo dataframe and include the 5end and 3end annotations 
#         if igRep.end5:     
#             if '5end' not in cloneAnnot.columns:
#                 cloneAnnot.loc[:, '5end'] = Series(valid5End, index=valid5End.keys())
#                 cloneAnnot.loc[:, '5endPrimer'] = Series(primer5End, index=primer5End.keys())
#                 cloneAnnot.loc[:, '5endIndel'] = Series(indel5End, index = indel5End.keys())
# 
#         if igRep.end3:    
#             if '3end' not in cloneAnnot.columns:        
#                 cloneAnnot.loc[:, '3end'] = Series(valid3End, index=valid3End.keys())
#                 cloneAnnot.loc[:, '3endPrimer'] = Series(primer3End, index=primer3End.keys())
#                 cloneAnnot.loc[:, '3endIndel'] = Series(indel3End, index = indel3End.keys())
# 
#         ## add columns of the stop codon location 
#         if 'stopat5end' not in cloneAnnot.columns:
#             df1 = DataFrame.from_dict(stopCodonPos, orient='index')
#             df1.columns = ['stopat5end', 'stopat3end', 'stopinmiddle']
#             cloneAnnot = concat([cloneAnnot, df1], axis=1)
# 
#         ## print statsitics and final processed data 
#         if (len(fr1NotAtBegin) > 0):
#             print("%d clones have FR1 start not equal to query start (Excluded)" % (len(fr1NotAtBegin)))
#             examples = random.choice(range(len(fr1NotAtBegin)), 10)
#             for i in examples:
#                 print(fr1NotAtBegin[i])
#             fr1NotAtBegin = None 
#         if (len(endsWithStopCodon) > 0):
#             print("%d clones contain a stop codon " % (len(endsWithStopCodon)))
#             examples = random.choice(range(len(endsWithStopCodon)), 10)
#             for i in examples:
#                 print(endsWithStopCodon[i])
#             endsWithStopCodon = None   
#         if (updatedStopCodon > 0):
#             print("The stopcodon flag was updated for %d sequencs " % (updatedStopCodon))          
#         if (len(fr4NotAsExpected) > 0):
#             print("%d clones do not have an expected FR4 clones (%s) " % (len(fr4NotAsExpected), FR4_CONSENSUS))
#             examples = random.choice(range(len(fr4NotAsExpected)), 10)
#             for i in examples:
#                 print(fr4NotAsExpected[i])
#             fr4NotAsExpected = None
#         if (len(noFR4) > 0):
#             print("%d clones do not have FR4 " % (len(noFR4)))
#             examples = random.choice(range(len(noFR4)), 10)
#             for i in examples:
#                 print(noFR4[i])
#             noFR4 = None
#             
#         igRep.cloneSeqs = DataFrame(transSeqs, columns=['queryid', 'germline', 'fr1', 'cdr1', 'fr2', 'cdr2',
#                                                   'fr3', 'cdr3', 'fr4'])
#         igRep.cloneSeqs.index = igRep.cloneSeqs.queryid
#         del igRep.cloneSeqs['queryid']
#         igRep.cloneSeqs.to_csv(igRep.cloneSeqFile, sep='\t', header=True, index=True)
#         
#         # export the CDR/FR annotation to a file
#         cloneAnnot.to_csv(igRep.cloneAnnotFile, sep='\t', header=True, index=True)    
#         print("Text file has been written to " + igRep.cloneAnnotFile) 
#         sys.stdout.flush()
#         gc.collect()
#         
#         
        
        
        
        