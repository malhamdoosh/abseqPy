import sys

from multiprocessing import Process

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from abseqPy.IgRepertoire.igRepUtils import findBestMatchedPattern, calMaxIUPACAlignScores
from abseqPy.logger import printto, LEVEL


class PrimerWorker(Process):
    def __init__(self, procCounter, fr4cut, trim5end,
                 trim3end, actualQstart, end5, end3,
                 end5offset, tasks, exitQueue, resultsQueue, stream=None):
        super(PrimerWorker, self).__init__()
        self.procCounter = procCounter
        self.fr4cut = fr4cut
        self.trim5end = trim5end
        self.trim3end = trim3end
        self.actualQstart = actualQstart
        self.end5 = end5
        self.end3 = end3
        self.end5offset = end5offset
        self.taskQueue = tasks
        self.exitQueue = exitQueue
        self.resultsQueue = resultsQueue
        self.firstJobTaken = False
        self.maxPrimer5Length, self.primer5sequences = _parsePrimerFile(self.end5, stream=stream)
        self.maxPrimer3Length, self.primer3sequences = _parsePrimerFile(self.end3, stream=stream)
        self.stream = stream

    def run(self):
        while True:
            nextTask = self.taskQueue.get()
            if nextTask is None:
                printto(self.stream, self.name + " process has stopped.")
                self.exitQueue.put("exit")
                break

            try:
                recs = []
                if not self.firstJobTaken:
                    printto(self.stream, self.name + " process commenced a new task ... ")
                    self.firstJobTaken = True
                for record, qsRec in zip(nextTask[0], nextTask[1]):
                    qsRec['queryid'] = record.id
                    recs.append(_matchClosestPrimer(qsRec, record, self.actualQstart, self.trim5end,
                                                    self.trim3end, self.end5offset, self.fr4cut, self.maxPrimer5Length,
                                                    self.maxPrimer3Length, self.primer5sequences,
                                                    self.primer3sequences))
                self.resultsQueue.put(recs)
                self.procCounter.increment(len(recs))
            except Exception as e:
                printto(self.stream, "An error as occurred while processing " + self.name + " with error {}".format(
                    str(e)
                ), LEVEL.EXCEPT)
                self.resultsQueue.put(None)
                continue
        return


def _matchClosestPrimer(qsRec, record, actualQstart, trim5end, trim3end, end5offset, fr4cut,
                        maxPrimer5Length, maxPrimer3Length, primer5seqs, primer3seqs):
    if qsRec['strand'] == 'reversed':
        record = SeqRecord(record.seq.reverse_complement(), id=record.id, name="", description="")

    vh = record.seq

    unexpected5 = unexpected3 = 0

    if primer5seqs:
        primer = str(vh[max(0, end5offset):max(0, end5offset) + maxPrimer5Length])
        try:
            qsRec['5endPrimer'], qsRec['5endMismatchIndex'], qsRec['5endIndelIndex'], _, _ = \
                findBestMatchedPattern(primer, primer5seqs)
        except Exception as e:
            # print("ARGH: something went wrong!" + str(e.message))
            unexpected5 += 1
            pass

    if primer3seqs:
        primer = str(vh[-1 * maxPrimer3Length:])
        try:
            qsRec['3endPrimer'], qsRec['3endMismatchIndex'], qsRec['3endIndelIndex'], _, _ = \
                findBestMatchedPattern(primer, primer3seqs)
        except Exception as e:
            # print("DEBUG: something went wrong! {}".format(str(e.message)))
            unexpected3 += 1
            pass

    # finish
    return qsRec, unexpected5, unexpected3


def _parsePrimerFile(primerFile, stream=None):
    if primerFile:
        primerids = []
        primerLengths = []
        primerSequences = []
        for rec in SeqIO.parse(primerFile, "fasta"):
            primerLengths.append(len(rec.seq))
            primerids.append(rec.id)
            primerSequences.append(str(rec.seq).upper())

        maxScores = calMaxIUPACAlignScores(primerSequences)

        if len(set(primerLengths)) != 1:
            printto(stream, "WARNING: Provided primer file {} has primers with different length. "
                            "Analysis assumes uniform primer length"
                  .format(primerFile), LEVEL.WARN)
        return max(primerLengths), zip(primerids, primerSequences, maxScores)

    return None, None
