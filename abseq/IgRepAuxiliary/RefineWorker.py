'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import sys
import os

from multiprocessing import Process
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from numpy import isnan, nan

from abseq.config import FR4_CONSENSUS, FR4_CONSENSUS_DNA
from abseq.IgRepertoire.igRepUtils import extractProteinFrag, \
    findBestAlignment, extractCDRsandFRsProtein, extractCDRsandFRsDNA, calMaxIUPACAlignScores, findBestMatchedPattern
from abseq.IgRepAuxiliary.IgBlastWorker import convertCloneRecordToOrderedList


class RefineWorker(Process):
    def __init__(self, procCounter, igdb, chain, actualQstart,
                 fr4cut, trim5End, trim3End, refineFlagNames):
        super(RefineWorker, self).__init__()
        self.procCounter = procCounter
        # we could have done the J gene length map reading earlier, but this encapsulates it better.
        # you'll only be reading N number of times where N = num.threads
        self.jGeneLengthMap = _generateJMap(igdb, chain)
        self.chain = chain
        self.actualQstart = actualQstart
        self.fr4cut = fr4cut
        self.trim5End = trim5End
        self.trim3End = trim3End if type(trim3End) == int else _parse3EndSeqs(trim3End)
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
                                                  self.chain, self.jGeneLengthMap, flags)
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
                          trim5End, trim3End, chain, jGeneLengthMap, flags):
    seqs = [record.id, qsRec['vgene']]

    try:
        if qsRec['strand'] == "reversed":
            record = SeqRecord(record.seq.reverse_complement(), id=record.id,
                               name="", description="")

        record = record[trim5End:]
        # if trim3End was a user provided int, use it to cut the sequence, or else it will be a list of seqs
        if type(trim3End) == int:
            record = record[:(len(record) - trim3End)]

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

            fr4start, fr4end, gapped = findBestAlignment(searchRegion, FR4_CONSENSUS[chain])

            if not gapped and fr4start != -1 and fr4end != -1 and fr4end > fr4start:
                # FR4 start
                qsRec['fr4.start'] = (fr4start - 1) * 3 + qsRec['fr3.end'] + 1

                # CDR3                
                qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
            else:
                # try to use the DNA consensus
                searchRegion = str(record.seq)[int(qsRec['fr3.end']):]
                fr4start, fr4end, gapped = findBestAlignment(searchRegion, FR4_CONSENSUS_DNA[chain], True)
#                 print ("DNA", searchRegion, fr4start, fr4end)  
                if fr4start != -1 and fr4end != -1 and fr4end > fr4start:
                    # FR4 start
                    qsRec['fr4.start'] = qsRec['fr3.end'] + fr4start

                    # CDR3                
                    qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
                    flags['CDR3dna'] += [record.id]
                else:
                    # TODO: check this case
                    qsRec['cdr3.end'] = qsRec['jqend']
                    # qsRec['fr4.end'] = len(record.seq)
                    # qsRec['fr4.start'] =  len(record.seq)

            # FR4 end
            # Check whether to cut the Ig sequence after FR4 or not
            if not fr4cut:
                # if trim3End wasn't of type int, it was a file with seqs to match
                if type(trim3End) != int:
                    _, _, _, relativeFR4EndPosition, _ = \
                        findBestMatchedPattern(str(record.seq)[int(qsRec['cdr3.end']):], trim3End, extend5end=True)
                    if relativeFR4EndPosition == -1:
                        # bad alignment - flag it and fall back to Jend
                        flags['FR4endless'] += [record.id]
                        qsRec['fr4.end'] = _predictFR4end(jGeneLengthMap, qsRec, flags, record)
                    else:
                        qsRec['fr4.end'] = qsRec['cdr3.end'] + relativeFR4EndPosition
                        if qsRec['fr4.end'] < qsRec['jqend']:
                            # if FR4 is somehow misaligned, flag it and fall back to Jend
                            # this may happen if: for example, the provided sequence happen to also appear
                            # before J germline ends, and it matches that first.
                            flags['FR4cutEarly'] += [record.id]
                            qsRec['fr4.end'] = _predictFR4end(jGeneLengthMap, qsRec, flags, record)
                else:
                    qsRec['fr4.end'] = len(record.seq)
            else:
                qsRec['fr4.end'] = _predictFR4end(jGeneLengthMap, qsRec, flags, record)

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


def _generateJMap(database, chain, species='human'):
    lengthMap = {}
    pathToDB = os.path.join(os.path.abspath(os.path.expandvars(database)), 'imgt_{}_ig{}j'.format(species, chain[0]))
    for rec in SeqIO.parse(pathToDB, 'fasta'):
        lengthMap[rec.id.strip()] = len(rec.seq.strip())
    return lengthMap


def _parse3EndSeqs(seqs):
    """
    transform list of seqs to expected format by findBestMatchedPattern
    :param seqs: raw string sequences
    :return: zippped(ids, seqs, maxUIPACScores)
    """
    targetids = [str(i) for i in range(len(seqs))]
    maxScores = calMaxIUPACAlignScores(seqs)
    return zip(targetids, seqs, maxScores)


def _predictFR4end(jGeneLengthMap, qsRec, flags, record):
    """
    predict FR4end by extending J germline end index
    :param jGeneLengthMap: map of J gene name to its length
    :param qsRec: cloneAnnot row in type dict
    :param flags: refinement flags
    :param record: original sequence (raw)
    :return: predicted FR4.end in nt position 1-based index, else np.nan if cannot predict.
    """
    # FR4end = num of remaining unmapped J nt + queryJEnd
    #        = (JGeneLength - subjectJEnd) + queryJEnd
    #        = (JGeneLength - (subjectJStart + (number of mapped J nt) - 1) ) + queryJEnd
    #        = (JGeneLength - (subjectJStart + (queryJEnd - queryJStart + 1) - 1) ) + queryJEnd
    #        = (JGeneLength - subjectJStart -(queryJEnd - queryJStart + 1) + 1) + queryJEnd
    #        = (JGeneLength - subjectJStart - queryJEnd + queryJStart - 1 + 1) + queryJEnd
    #        = JGeneLength - subjectJStart - queryJEnd + queryJStart + queryJEnd
    #        = JGeneLength - subjectJStart + queryJStart
    # why FR4end = num unmapped J nt + queryJend? - igblast rather cut J gene instead of allowing
    # mismatches => no "extend 3' end" option. So we add remaining unmapped j nt to jqueryend to get
    # end of FR4 (end of whole J gene)
    if type(qsRec['jgene']) == str and qsRec['jgene'] != '' and qsRec['jgene'].strip() in jGeneLengthMap:
        # if jgene is not nan, we can assume that jstart(subject), jqend, jqstart are all defined
        predictedFR4End = jGeneLengthMap[qsRec['jgene'].strip()] - qsRec['jstart'] + qsRec['jqstart']
        if predictedFR4End < 0 or predictedFR4End < qsRec['jqend']:
            predictedFR4End = qsRec['jqend']
            # this usually doesn't happen
            flags['FR4PredictedError'] += [record.id]
        if predictedFR4End > len(record.seq):
            # this will happen if your sequences end before J gene ends
            flags['FR4AbruptEnd'] += [record.id]
        return min(len(record.seq), predictedFR4End)
    return nan
