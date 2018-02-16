'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

from multiprocessing import Process
import sys
from Bio.SeqRecord import SeqRecord
from config import FR4_CONSENSUS, FR4_CONSENSUS_DNA
from IgRepertoire.igRepUtils import extractProteinFrag,\
    findBestAlignment, extractCDRsandFRsProtein, extractCDRsandFRsDNA
from numpy import isnan
import traceback
from IgBlastWorker import convertCloneRecordToOrderedList


class RefineWorker(Process):
    def __init__(self, procCounter, chain, actualQstart, 
                 fr4cut, trim5End, trim3End, refineFlagNames):
        super(RefineWorker, self).__init__()
        self.procCounter = procCounter
        self.chain = chain
        self.actualQstart = actualQstart
        self.fr4cut = fr4cut
        self.trim5End = trim5End
        self.trim3End = trim3End
        self.refineFlagNames = refineFlagNames
        self.tasksQueue = None
        self.exitQueue = None
        self.resultsQueue = None   
        self.firstJobTaken = False 
    def run(self):
        print(self.name + " process is now ready to start a new job ...")
        while True:            
            nextTask = self.tasksQueue.get()
            # poison pill check            
            if (nextTask is None):
                print(self.name + " process has stopped." )
                self.exitQueue.put("exit")
                break
            try:
                if not self.firstJobTaken: 
                    print(self.name + " process commenced a new task ... ")
                    self.firstJobTaken = True
                qsRecs = []
                seqsAll = []
                flags = {} 
                for f in self.refineFlagNames:
                    flags[f] = []
                for (record, qsRec) in zip(nextTask[0], nextTask[1]):  
                    seqs = refineCloneAnnotation(qsRec, record, 
                                                  self.actualQstart, self.fr4cut, 
                                                  self.trim5End, self.trim3End,
                                                  self.chain, flags) 
                    # out-of-frame clones are excluded
                    if qsRec['v-jframe'] != 'Out-of-frame':
                        refineInFramePrediction(qsRec, record, self.actualQstart, flags)
                    
                    # append the FR and CDR protein clones
                    qsRec['queryid'] = record.id
                    qsRecs.append(convertCloneRecordToOrderedList(qsRec, self.chain))
                    seqsAll.append(seqs)        
                self.procCounter.increment(len(qsRecs))                         
                self.resultsQueue.put((qsRecs, seqsAll, flags))
#                 print('\t%d clones have been processed ... ' % (len(nextTask)))
                sys.stdout.flush()        
            except Exception as e:
                print("An error occurred while processing " + self.name)
                print(e)
                self.resultsQueue.put(None)
                continue                       
#             print("process has completed a run... " + self.name) 
        return


def refineCloneAnnotation(qsRec, record, actualQstart, fr4cut, 
                          trim5End, trim3End, chain, flags):
    seqs = [record.id, qsRec['vgene']]

    try:
        if qsRec['strand'] == "reversed":
            record = SeqRecord(record.seq.reverse_complement(), id=record.id,
                               name="", description="")
        record = record[trim5End:(len(record) - trim3End)]

        # grab the beginning of the VH clone
        if actualQstart > -1:
            # if user specified an actualQstart, use it. (parse args already converted it into 0-based)
            offset = actualQstart # zero-based
        else:
            # else, we find the offset by subtracting V query start with IgBLAST's v germline start position (1-index)
            offset = int(qsRec['vqstart'] - qsRec['vstart'])  # zero-based
        if offset < 0:
            offset = 0                   
        vh = record.seq[offset:]

        # check whether the VH sequence can be translated successfully
        if len(vh) % 3 != 0:
            vh = vh[:-1 * (len(vh) % 3)]
        protein = str(vh.translate())

        # check whether the start of the V gene is the same as the start of FR1
        if qsRec['vqstart'] != qsRec['fr1.start']:
            flags['fr1NotAtBegin'] += [record.id]
        qsRec['fr1.start'] = offset + 1

        # Identification of FR4 so that CDR3 can be defined
        if isnan(qsRec['fr4.end']):            
            searchRegion = extractProteinFrag(protein, qsRec['fr3.end'] + 1, -1, offset, trimAtStop=False)
            if searchRegion is None:
                raise Exception("ERROR: undefined search region to find FR3 consensus.")
            qsRec['cdr3.start'] = qsRec['fr3.end'] + 1
            fr4start, fr4end, gapped = findBestAlignment(searchRegion, 
                                                         FR4_CONSENSUS[chain])# , show=True
#             print ("Protein", searchRegion, fr4start, fr4end)                 
            if not gapped and fr4start != -1 and fr4end != -1 and fr4end > fr4start:
                qsRec['fr4.start'] = (fr4start - 1) * 3 + qsRec['fr3.end'] + 1
                # CDR3                
                qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
                # Check whether to cut the Ig sequence after FR4 or not
                if not fr4cut: 
                    qsRec['fr4.end'] = len(record.seq)  # fr4end * 3 + qsRec['fr3.end']
                else:
                    qsRec['fr4.end'] = qsRec['fr3.end'] + fr4end * 3
            else:                
                # try to use the DNA consensus
                searchRegion = str(record.seq)[int(qsRec['fr3.end']):]
                fr4start, fr4end, gapped = findBestAlignment(searchRegion, 
                                                             FR4_CONSENSUS_DNA[chain],
                                                             True)  # , show=True
#                 print ("DNA", searchRegion, fr4start, fr4end)  
                if fr4start != -1 and fr4end != -1 and fr4end > fr4start:
                    qsRec['fr4.start'] = qsRec['fr3.end'] + fr4start
                    # CDR3                
                    qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
                    if not fr4cut: 
                        qsRec['fr4.end'] = len(record.seq) 
                    else:
                        qsRec['fr4.end'] = qsRec['fr3.end'] + fr4end
                    flags['CDR3dna'] += [record.id]
                else:
                    # TODO: check this case
                    qsRec['cdr3.end'] = qsRec['jqend']
#                     qsRec['fr4.end'] = len(record.seq) 
#                     qsRec['fr4.start'] =  len(record.seq)
        # Extract the CDR and FR protein sequences
        (protein, tmp) = extractCDRsandFRsProtein(protein, qsRec, offset)
        seqs += tmp
        if seqs[-1][:4] != FR4_CONSENSUS[chain][:4]:
            flags['fr4NotAsExpected'] += [record.id]
        if seqs[-1] == '':
            flags['noFR4'] += [record.id]

        # Extract the CDR and FR nucleotide sequences
        # COMMENTED OUT ON:  Fri Feb 16 10:36:41 AEDT 2018 BY JIAHONG - REASON: tmp var not used
        # is also generating "clones not partitioned correctly although the protein seqs are"
        #tmp = extractCDRsandFRsDNA(str(record.seq), qsRec)
        # TODO: consider adding the DNA sequences
        if '*' in protein:
            flags['endsWithStopCodon'] += [record.id]                     
            # update the StopCodon value if it was set to No
            if qsRec['stopcodon'] == 'No':
                flags['updatedStopCodon'] += [record.id]
                qsRec['stopcodon'] = 'Yes' 
        
        # update the annotation fields with the new calculated values                        
        gaps = abs(qsRec['vqstart'] - qsRec['vstart']) - offset                
        mismatches = qsRec['vstart'] - 1
        if qsRec['vstart'] > qsRec['vqstart'] and gaps > 0:
            mismatches -= gaps

        # Only update gaps if the actual query start position is known 
        if gaps > 0:     
            qsRec['fr1.gaps'] += gaps                    
            qsRec['vgaps'] += gaps

        # if igblast ignores mismatches at the beginning ==> update
        if mismatches > 0:
            qsRec['fr1.mismatches'] += mismatches
            qsRec['vmismatches'] += mismatches
            qsRec['vstart'] -= mismatches
            qsRec['vqstart'] -= mismatches
        # TODO: update gaps and mismatches in FR4 and CDR3 based on D and J germlines
        # TODO: update the start and end fields based on the trim5End
    except Exception as e:    
        if "partitioning" in str(e):
            flags['partitioning'] += [record.id]            
#         print("ERROR: exception in the Clone Annotation Refinement method")
#         print(record.id)
#         print(str(vh))
#         print(qsRec)
#         print(protein)        
#         traceback.print_exc(file=sys.stdout)
#         raise e
    return seqs


def refineInFramePrediction(qsRec, record, actualQstart, flags):
        inframe = True    
        # check the the v-jframe value is not NA       
        if (qsRec['v-jframe'] == 'N/A' or
            (type(qsRec['v-jframe']) != type('str') and isnan(qsRec['v-jframe']))):
            flags['updatedInFrameNA'] += [record.id]
            inframe = False
        # the query clone is not in concordance with the start of the germline gene
        offset = qsRec['vqstart'] - qsRec['vstart'] + 1 # 1-based 
        if (inframe and (offset < 1 or 
            (actualQstart != -1 and (offset - 1 - actualQstart) % 3 != 0))):            
            inframe = False
            flags['updatedInFrameConc'] += [record.id]           
        # if no CDR3 or FR4 ==> Out-of-frame
        if (inframe and (isnan(qsRec['fr4.start']) or 
            isnan(qsRec['fr4.end'])  or isnan(qsRec['cdr3.start']) or 
           qsRec['cdr3.start'] >= qsRec['cdr3.end'])):
            inframe = False
            flags['updatedInFrameNo3or4'] += [record.id]                        
        # doesn't start/end properly .. not multiple of 3
        if (inframe and ((qsRec['fr4.end'] - qsRec['fr1.start'] + 1) % 3 != 0 or 
            (actualQstart != -1 and 
             ((qsRec['fr4.end'] - actualQstart) % 3 != 0)))):
            inframe = False
            flags['updatedInFrame3x'] += [record.id]
#             print(str(record.seq))
        # indels (gaps) in FRs or CDRs cause frame-shift ==> out-of-frame
        if (inframe and ( 
            qsRec['fr1.gaps'] % 3 != 0 or 
            qsRec['fr2.gaps'] % 3 != 0 or
            qsRec['fr3g.gaps'] % 3 != 0 or
            qsRec['cdr1.gaps'] % 3 != 0 or 
            qsRec['cdr2.gaps'] % 3 != 0 or
            (not isnan(qsRec['cdr3g.gaps']) and qsRec['cdr3g.gaps'] % 3 != 0)
            )):
            inframe = False            
            flags['updatedInFrameIndel'] += [record.id]
        # FR1 start is not aligned with query start 
#             if (not inframe or qsRec['vqstart'] != qsRec['fr1.start']):
#                 inframe = False
        if not inframe:
            qsRec['v-jframe'] =  'Out-of-frame'
            flags['updatedInFrame'] += [record.id]
        


