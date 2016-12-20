from Bio import SeqIO
from pandas.core.frame import DataFrame
from numpy import  random
from multiprocessing import Queue, Process, Manager, Value, Lock
from IgRepAuxiliary.RefineWorker import RefineWorker
import sys
from math import ceil
from config import MEM_GB
import gc
from IgRepAuxiliary.IgBlastWorker import getAnnotationFields

def loadRefineFlagInfo():
    refineFlagNames = ['fr1NotAtBegin', 'endsWithStopCodon', 
             'fr4NotAsExpected', 'updatedStopCodon', 'noFR4', 'CDR3dna',
             'partitioning',
             'updatedInFrame' , 'updatedInFrameNA', 'updatedInFrameConc',
             'updatedInFrameNo3or4', 'updatedInFrame3x' ,
             'updatedInFrameIndel']
    refineFlagMsgs = {}
    refineFlagMsgs['fr1NotAtBegin'] = "{:,} clones have FR1 start not equal to query start (Excluded)" 
    refineFlagMsgs['endsWithStopCodon'] = "{:,} clones contain a stop codon "
    refineFlagMsgs['updatedStopCodon'] = "The stopcodon flag was updated for {:,} clones "
    refineFlagMsgs['fr4NotAsExpected'] = "The FR4 of {:,} clones do not start as expected"
    refineFlagMsgs['noFR4'] = "{:,} clones do not have FR4 "
    refineFlagMsgs['updatedInFrame'] = "The v-j frame rearrangement status has been corrected for {:,} clones "
    refineFlagMsgs['updatedInFrameNA'] = "{:,} clones have undefined in-frame status"
    refineFlagMsgs['updatedInFrameConc'] = "{:,} clones show discordance between the query and v gene starts"
    refineFlagMsgs['updatedInFrameNo3or4'] = "{:,} clones have no CDR3 or FR4"
    refineFlagMsgs['updatedInFrame3x'] = "{:,} clones are not multiple of 3 "
    refineFlagMsgs['updatedInFrameIndel'] = "{:,} clones have indels in one of the FRs or CDRs"
    refineFlagMsgs['CDR3dna'] = "The CDR3 of {:,} clones was determined using DNA consensus"
    refineFlagMsgs['partitioning'] = "{:,} clones were partitioned incorrectly."
    return (refineFlagNames, refineFlagMsgs)

def refineClonesAnnotation(outDir, sampleName, cloneAnnotOriginal, readFile, format, 
                            actualQstart, chain, fr4cut, 
                            trim5End, trim3End,
                            seqsPerFile, threads):
        print("Clone annotation and in-frame prediction are being refined ...")        
        seqsPerFile = 100
        cloneAnnot = cloneAnnotOriginal.copy()        
        queryIds = cloneAnnot.index#[4200000:]
        (refineFlagNames, refineFlagMsgs) = loadRefineFlagInfo()
#         manager = Manager()        
        try:    
            # process clones from the FASTA/FASTQ file
            records = SeqIO.index(readFile, format)    
            print("\t " +  format + " index created and refinement started ...")    
            ### Parallel implementation of the refinement
            noSeqs = len(queryIds)
            totalTasks = int(ceil(noSeqs * 1.0 / seqsPerFile)) 
            tasks = Queue()      
            exitQueue = Queue()
            resultsQueue = Queue()        
            procCounter = ProcCounter(noSeqs)    
            if (threads > totalTasks):
                threads = totalTasks     
            if MEM_GB < 16:
                threads = 2  
            # Initialize workers 
            workers = []        
            for i in range(threads):
                w = RefineWorker(procCounter, chain, actualQstart, 
                                 fr4cut, trim5End, trim3End, refineFlagNames)
                w.tasksQueue = tasks
                w.exitQueue = exitQueue  
                w.resultsQueue = resultsQueue    
                workers.append(w)
                w.start()       
                sys.stdout.flush()          
            # adding jobs to the tasks queue with subsets of query IDs
            if (totalTasks > 1): 
                for i in range(totalTasks):
                    ids = queryIds[i*seqsPerFile:(i+1)*seqsPerFile]
#                     print(ids)
#                     print(records[0])
                    recs = map(lambda x: records[x], ids)
                    qsRecs = map(lambda x: cloneAnnot.loc[x].to_dict(), ids)
                    tasks.put((recs, qsRecs))                
            else:                
                recs = map(lambda x: records[x], queryIds)
                qsRecs = map(lambda x: cloneAnnot.loc[x].to_dict(), queryIds)
                tasks.put((recs, qsRecs))
            # Add a poison pill for each worker
            for i in range(threads + 10):
                tasks.put(None)       
            # Wait all process workers to terminate                
            i = 0 
            while i < threads:    
                m = exitQueue.get()
                if m == "exit":
                    i += 1
            print("All workers have completed their tasks successfully.")
            # Collect results
            print("Results are being collated from all workers ...")  
            # invoking the result collection method   
            cloneAnnotList, transSeqs, flags = collectRefineResults(resultsQueue, totalTasks, 
                       noSeqs, refineFlagNames)
            ### End of parallel implementation                  
            sys.stdout.flush()           

            print("\tResults were collated successfully.")
            # print refine flags 
            printRefineFlags(flags, records, refineFlagNames, refineFlagMsgs)
            writeRefineFlags(flags, records, refineFlagNames, refineFlagMsgs, 
                             outDir, sampleName)
#             sys.exit()
        except Exception as e:
            print("Something went wrong during the refinement process!")
            raise e
        finally:
            for w in workers:
                w.terminate()     
            records.close()
        # Create new data frame of clone annotation
        cloneAnnot = DataFrame(cloneAnnotList, columns=getAnnotationFields(chain))
        cloneAnnot.index = cloneAnnot.queryid
        del cloneAnnot['queryid']
        gc.collect()        
        # Create data frame of FR and CDR sequences
        cols = ['queryid', 'germline', 'fr1', 'cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3', 'fr4']
        cloneSeqs = DataFrame(transSeqs, columns=cols)
        for col in cols:
            cloneSeqs.loc[:, col] = cloneSeqs[col].map(str)
        cloneSeqs.index = cloneSeqs.queryid
        del cloneSeqs['queryid']        
        return (cloneAnnot, cloneSeqs)
        
def collectRefineResults(resultsQueue, totalTasks, 
                       noSeqs, refineFlagNames):    
    total = 0
    cloneAnnot = []
    transSeqs = []
    flags = {}
    for f in refineFlagNames:
        flags[f] = []
    while totalTasks:                
        result = resultsQueue.get()
        totalTasks -= 1                           
        if (result is None):
            continue        
#         print(total, resultsQueue.qsize())      
        qsRecsOrdered, seqs, flagsi = result        
        # update relevant annotation fields
        cloneAnnot += qsRecsOrdered
        transSeqs +=  seqs  
        # update flags 
        for f in refineFlagNames:
            flags[f] += flagsi[f]
        total += len(qsRecsOrdered)
#         if (len(qsRecsOrdered) < 100):
#             print(len(qsRecsOrdered))
        if (total % 50000 == 0):
            print('\t%d/%d records have been collected ... ' % (total, noSeqs))
            sys.stdout.flush()
    print('\t%d/%d records have been collected ... ' % (total, noSeqs))
    return cloneAnnot, transSeqs, flags        
          
def printRefineFlags(flags, records, refineFlagNames, refineFlagMsgs):
    ## print statistics and a few of the flagged clones 
    for f in refineFlagNames:
        if (len(flags[f]) > 0):
            print(refineFlagMsgs[f].format(len(flags[f])))
            examples = random.choice(range(len(flags[f])), min(3, len(flags[f])))
            for i in examples:
                print(">" + flags[f][i])
                print(str(records[flags[f][i]].seq))
            print      
             
def writeRefineFlags(flags, records, refineFlagNames, 
                     refineFlagMsgs, outDir, sampleName):
    with open(outDir + sampleName + "_refinement_flagged.txt", 'w') as out:
        for f in refineFlagNames:
            if (len(flags[f]) > 0):
                out.write("# " + refineFlagMsgs[f].format(len(flags[f])) + "\n")                
                for i in range(len(flags[f])):
                    out.write(">" + flags[f][i] + "\n")
                    out.write(str(records[flags[f][i]].seq) + "\n")
                out.write("\n")  
     
     
class ProcCounter(object):
    def __init__(self, noSeqs, initval=0):
        self.val = Value('i', initval)
        self.noSeqs = noSeqs
        self.lock = Lock()

    def increment(self, val=1):
        with self.lock:
            self.val.value += val
            if (self.val.value % 50000 == 0 or self.noSeqs == self.val.value):
                print('\t%d/%d records have been processed ... ' % (self.val.value, self.noSeqs))

    def value(self):
        with self.lock:
            return self.val.value
        
        
           
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
# #                 if offset <= igRep.actualQstart:
# #                     continue
#                 if (qsRec['strand'] == "reversed"):
#                     record.seq = record.seq.reverse_complement()
#                 # grab the beginning of the VH clone  
#                 if igRep.actualQstart > -1:
#                     offset = igRep.actualQstart # zero-based
#                 else:                                  
#                     offset = int(qsRec['vqstart'] - qsRec['vstart'])  # zero-based
#                 if  offset < 0:
#                     offset = 0                   
#                 vh = record.seq[offset:]                       
#                 # check whether the VH clone can be translated successfully
#                 if len(vh) % 3 != 0:
#                     vh = vh[:-1 * (len(vh) % 3)]#                   
#                 protein = str(vh.translate())
#                 
#                 if qsRec['vqstart'] != qsRec['fr1.start']:
#                     fr1NotAtBegin += [record.id]
# 
#                 # FR1             
#                 seqs.append(extractProteinFrag(protein, offset,
#                                                qsRec['fr1.end'], offset))
#                 # CDR1
#                 seqs.append(extractProteinFrag(protein, qsRec['cdr1.start'],
#                                                qsRec['cdr1.end'], offset))
#                 # FR2
#                 seqs.append(extractProteinFrag(protein, qsRec['fr2.start'],
#                                                qsRec['fr2.end'], offset))
#                 # CDR2
#                 seqs.append(extractProteinFrag(protein, qsRec['cdr2.start'],
#                                                qsRec['cdr2.end'], offset))
#                 # FR3
#                 seqs.append(extractProteinFrag(protein, qsRec['fr3.start'],
#                                                qsRec['fr3.end'], offset))
#                 # Identification of FR4 so that CDR3 can be defined 
#                 if isnan(qsRec['fr4.end']):
#                     fr4start, fr4end = findBestAlignment(
#                                 extractProteinFrag(protein, qsRec['fr3.end'] + 1,
#                              - 1, offset, trimAtStop=False), FR4_CONSENSUS)#                     
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
#                                                    qsRec['cdr3.end'], offset))
#                 seqs.append(extractProteinFrag(protein, qsRec['fr4.start'],
#                                                qsRec['fr4.end'], offset))
#                 ## Check whether to cut the Ig clone after FR4 or not
#                 if igRep.fr4cut:
#                     try:                        
#                         protein = ''.join(seqs[2:])
#                     except:
#                         pass
#                     try:
#                         vh = record.seq[offset:int(qsRec['fr4.end'])]
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
#                         if offset + igRep.end5offset >= 0:
#                             prim = str(record.seq[offset + igRep.end5offset: offset + L5 + igRep.end5offset])
#                         else:
#                             prim = str(record.seq[: offset + L5 + igRep.end5offset])
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
#                 cloneAnnot.set_value(record.id, 'fr1.start', offset+1)                
#                 gaps = abs(qsRec['vqstart'] - qsRec['vstart']) - offset                
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
        
        
        
        