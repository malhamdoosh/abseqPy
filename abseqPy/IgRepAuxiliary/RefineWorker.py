'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

from multiprocessing import Process
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from numpy import isnan, nan

from abseqPy.config import FR4_CONSENSUS, FR4_CONSENSUS_DNA
from abseqPy.IgRepertoire.igRepUtils import extractProteinFrag, \
    findBestAlignment, extractCDRsandFRsProtein, calMaxIUPACAlignScores, findBestMatchedPattern
from abseqPy.IgRepAuxiliary.IgBlastWorker import convertCloneRecordToOrderedList
from abseqPy.logger import LEVEL, printto


class RefineWorker(Process):
    def __init__(self, procCounter, chain, actualQstart,
                 fr4cut, trim5End, trim3End, refineFlagNames, stream=None):
        super(RefineWorker, self).__init__()
        self.procCounter = procCounter
        self.chain = chain
        self.actualQstart = actualQstart
        self.fr4cut = fr4cut
        self.trim5End = trim5End
        self.trim3End = trim3End if isinstance(trim3End, int) else _parse3EndSeqs(trim3End)
        self.refineFlagNames = refineFlagNames
        self.tasksQueue = None
        self.exitQueue = None
        self.resultsQueue = None
        self.firstJobTaken = False
        self.stream = stream

    def run(self):
        printto(self.stream, self.name + " process is now ready to start a new job ...")
        while True:
            nextTask = self.tasksQueue.get()
            # poison pill check            
            if nextTask is None:
                printto(self.stream, self.name + " process has stopped.")
                self.exitQueue.put("exit")
                break
            try:
                if not self.firstJobTaken:
                    printto(self.stream, self.name + " process commenced a new task ... ")
                    self.firstJobTaken = True
                qsRecs = []
                seqsAll = []
                recordLengths = defaultdict(_defaultdefaultInt)
                flags = {}
                for f in self.refineFlagNames:
                    flags[f] = []
                for (record, qsRec) in zip(nextTask[0], nextTask[1]):
                    seqs = refineCloneAnnotation(qsRec, record,
                                                 self.actualQstart, self.chain, self.fr4cut,
                                                 self.trim5End, self.trim3End, flags, stream=self.stream)
                    # out-of-frame clones are excluded
                    if qsRec['v-jframe'] != 'Out-of-frame':
                        stillInFrame = refineInFramePrediction(qsRec, record, self.actualQstart,
                                                               flags, stream=self.stream)
                        if stillInFrame:
                            _recordFRLength(qsRec, recordLengths)

                    # append the FR and CDR protein clones
                    qsRec['queryid'] = record.id
                    qsRecs.append(convertCloneRecordToOrderedList(qsRec, self.chain))
                    seqsAll.append(seqs)
                self.procCounter.increment(len(qsRecs))
                self.resultsQueue.put((qsRecs, seqsAll, flags, recordLengths))
            except Exception as e:
                printto(self.stream, "An error occurred while processing " + self.name, LEVEL.EXCEPT)
                self.resultsQueue.put(None)
                continue
        return


def refineCloneAnnotation(qsRec, record, actualQstart, chain, fr4cut,
                          trim5End, trim3End, flags, stream=None):
    seqs = [record.id, qsRec['vgene']]

    if qsRec['chain'] in ['VH', 'VK', 'VL']:
        chain = qsRec['chain']
    else:
        # convert hv, kv, lv and klv to VH, VK, or VL
        if chain == 'klv':
            if 'K' in qsRec['vgene']:
                chain = "VK"
            else:
                # we run a risk that it's not VL, but we honestly have no other clue to deduce what chain this is
                chain = "VL"
        else:
            # if user provided a chain type, just use it as-is, although we will need to convert it to VH, VK, or VL
            # instead of hv, kv, lv
            chain = chain[::-1].upper()

        printto(stream, "Chain had unknown type {}".format(qsRec['chain']), LEVEL.WARN)

    try:
        if qsRec['strand'] == "reversed":
            record = SeqRecord(record.seq.reverse_complement(), id=record.id, name="", description="")

        record = record[trim5End:]
        # if trim3End was a user provided int, use it to cut the sequence, or else it will be a list of seqs
        if isinstance(trim3End, int):
            record = record[:(len(record) - trim3End)]

        # grab the beginning of the VH clone
        if actualQstart > -1:
            # if user specified an actualQstart, use it. (parse args already converted it into 0-based)
            offset = actualQstart  # zero-based
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
            searchRegion = extractProteinFrag(protein, qsRec['fr3.end'] + 1, -1, offset, trimAtStop=False,
                                              stream=stream)
            if searchRegion is None:
                raise Exception("ERROR: undefined search region to find FR4 consensus.")

            qsRec['cdr3.start'] = qsRec['fr3.end'] + 1

            fr4start, fr4end, gapped = findBestAlignment(searchRegion, FR4_CONSENSUS[chain], dna=False)

            if not gapped and fr4start != -1 and fr4end != -1 and fr4end > fr4start:

                qsRec['fr4.start'] = (fr4start - 1) * 3 + qsRec['fr3.end'] + 1

                qsRec['cdr3.end'] = qsRec['fr4.start'] - 1

                # FR4end here is amino acid indexing, change to nt
                fr4end = qsRec['fr3.end'] + fr4end * 3

            else:
                # try to use the DNA consensus
                searchRegion = str(record.seq)[int(qsRec['fr3.end']):]
                fr4start, fr4end, gapped = findBestAlignment(searchRegion, FR4_CONSENSUS_DNA[chain], dna=True)

                if fr4start != -1 and fr4end != -1 and fr4end > fr4start:
                    qsRec['fr4.start'] = qsRec['fr3.end'] + fr4start

                    qsRec['cdr3.end'] = qsRec['fr4.start'] - 1
                    flags['CDR3dna'] += [record.id]

                    # FR4end is already in nt counts (don't need to * 3)
                    fr4end = qsRec['fr3.end'] + fr4end
                else:
                    # TODO: check this case
                    qsRec['cdr3.end'] = qsRec['jqend']

                    # Can't deduce FR4 !
                    fr4end = nan
                    # qsRec['fr4.end'] = len(record.seq)
                    # qsRec['fr4.start'] =  len(record.seq)

            # FR4 end
            # Check whether to cut the Ig sequence after FR4 or not
            if not fr4cut:
                # if trim3End wasn't of type int, it was a file with seqs to match
                if not isinstance(trim3End, int):
                    _, _, _, relativeFR4EndPosition, _ = \
                        findBestMatchedPattern(str(record.seq)[int(qsRec['cdr3.end']):], trim3End, extend5end=True)
                    if relativeFR4EndPosition == -1:
                        # bad alignment - flag it and fall back to consensus FR4end
                        flags['FR4endless'] += [record.id]
                        qsRec['fr4.end'] = fr4end
                    else:
                        qsRec['fr4.end'] = qsRec['cdr3.end'] + relativeFR4EndPosition
                        if qsRec['fr4.end'] < qsRec['jqend']:
                            # if FR4 is somehow misaligned, flag it and fall back to consensus FR4end
                            # this may happen if: for example, the provided sequence happen to also appear
                            # before J germline ends, and it matches that first.
                            flags['FR4cutEarly'] += [record.id]
                            qsRec['fr4.end'] = fr4end
                else:
                    # trim3End is int, use sliced string from earlier as FR4end
                    qsRec['fr4.end'] = len(record.seq)
            else:
                if isnan(fr4end):
                    qsRec['FR4PredictedError'] += [record.id]
                qsRec['fr4.end'] = fr4end

        # Extract the CDR and FR protein sequences
        (protein, tmp) = extractCDRsandFRsProtein(protein, qsRec, offset, stream=stream)
        seqs += tmp
        if seqs[-1][:4] != FR4_CONSENSUS[chain][:4]:
            flags['fr4NotAsExpected'] += [record.id]
        if seqs[-1] == '':
            flags['noFR4'] += [record.id]

        # Extract the CDR and FR nucleotide sequences
        # COMMENTED OUT ON:  Fri Feb 16 10:36:41 AEDT 2018 BY JIAHONG - REASON: tmp var not used
        # is also generating "clones not partitioned correctly although the protein seqs are"
        # tmp = extractCDRsandFRsDNA(str(record.seq), qsRec)
        # TODO: consider adding the DNA sequences
        if '*' in protein:
            flags['endsWithStopCodon'] += [record.id]
            # update the StopCodon value if it was set to No
            if qsRec['stopcodon'] == 'No':
                flags['updatedStopCodon'] += [record.id]
                qsRec['stopcodon'] = 'Yes'

                # update the annotation fields with the new calculated values
        gaps = int(abs(qsRec['vqstart'] - qsRec['vstart']) - offset)
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
    return seqs


def refineInFramePrediction(qsRec, record, actualQstart, flags, stream=None):
    inframe = True

    # check the the v-jframe value is not NA
    if (qsRec['v-jframe'] == 'N/A' or
            (not isinstance(qsRec['v-jframe'], str) and isnan(qsRec['v-jframe']))):
        flags['updatedInFrameNA'] += [record.id]
        inframe = False

    # the query clone is not in concordance with the start of the germline gene
    offset = int(qsRec['vqstart'] - qsRec['vstart']) + 1  # 1-based
    if (inframe and (offset < 1 or
                     (actualQstart != -1 and (offset - 1 - actualQstart) % 3 != 0))):
        inframe = False
        flags['updatedInFrameConc'] += [record.id]

    # if no CDR3 or FR4 ==> Out-of-frame
    if (inframe and (isnan(qsRec['fr4.start']) or
                     isnan(qsRec['fr4.end']) or isnan(qsRec['cdr3.start']) or
                     qsRec['cdr3.start'] >= qsRec['cdr3.end'])):
        inframe = False
        flags['updatedInFrameNo3or4'] += [record.id]

    # doesn't start/end properly .. not multiple of 3
    if (inframe and ((qsRec['fr4.end'] - qsRec['fr1.start'] + 1) % 3 != 0 or
                     (actualQstart != -1 and
                      ((qsRec['fr4.end'] - actualQstart) % 3 != 0)))):
        inframe = False
        flags['updatedInFrame3x'] += [record.id]
    # print(str(record.seq))

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
        qsRec['v-jframe'] = 'Out-of-frame'
        flags['updatedInFrame'] += [record.id]

    return inframe


def _parse3EndSeqs(seqs):
    """
    transform list of seqs to expected format by findBestMatchedPattern
    :param seqs: raw string sequences
    :return: zippped(ids, seqs, maxUIPACScores)
    """
    targetids = [str(i) for i in range(len(seqs))]
    maxScores = calMaxIUPACAlignScores(seqs)
    return zip(targetids, seqs, maxScores)


def _recordFRLength(qsRec, germlineConsensusLength):
    vgene = qsRec['vgene'].split('*')[0]
    jgene = qsRec['jgene'].split('*')[0]
    for region in ('fr1', 'fr2', 'fr3', 'fr4'):
        start = region + ".start"
        end = region + '.end'
        length = qsRec[end] - qsRec[start] + 1
        gene = vgene if region != 'fr4' else jgene
        if not isnan(length):
            germlineConsensusLength[gene][region][length] += 1
    return germlineConsensusLength


def _defaultInt():
    """
    to be pickle-able, this function cannot be lambda
    :return: equivalent to defaultdict(int)
    """
    return defaultdict(int)


def _defaultdefaultInt():
    """
    to be pickle-able, this function cannot be lambda
    :return: equivalent to defaultdict(lambda: defaultdict(int)
    """
    return defaultdict(_defaultInt)


