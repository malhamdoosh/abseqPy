import numpy as np

from multiprocessing import Process
from Bio.SeqRecord import SeqRecord

from abseq.IgRepertoire.igRepUtils import findBestMatchedPattern
from abseq.IgRepAuxiliary.primerAuxiliary import parsePrimerFile


class PrimerWorker(Process):
    def __init__(self, procCounter, fr4cut, trim5end,
                 trim3end, actualQstart, end5, end3,
                 end5offset, tasks, exitQueue, resultsQueue):
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
        self.maxPrimer5Length, self.primer5sequences = parsePrimerFile(self.end5)
        self.maxPrimer3Length, self.primer3sequences = parsePrimerFile(self.end3)

    def run(self):
        while True:
            nextTask = self.taskQueue.get()
            if nextTask is None:
                print(self.name + " process has stopped.")
                self.exitQueue.put("exit")
                break

            try:
                if not self.firstJobTaken:
                    print(self.name + " process commenced a new task ... ")
                    self.firstJobTaken = True
                for record, qsRec in zip(nextTask[0], nextTask[1]):
                    _matchClosestPrimer(qsRec, record, self.actualQstart, self.trim5end,
                                        self.trim3end, self.end5offset, self.fr4cut, self.maxPrimer5Length,
                                        self.maxPrimer3Length, self.primer5sequences, self.primer3sequences)
            except Exception as e:
                print("An error as occurred while processing " + self.name)
                print(e)
                self.resultsQueue.put(None)
                continue
        return


def _matchClosestPrimer(qsRec, record, actualQstart, trim5end, trim3end, end5offset, fr4cut,
                        maxPrimer5Length, maxPrimer3Length, primer5seqs, primer3seqs):
    if qsRec['strand'] == 'reversed':
        record = SeqRecord(record.seq.reverse_complement(), id=record.id, name="", description="")
    record = record[trim5end:(len(record) - trim3end)]

    if actualQstart > -1:
        # zero based (argparse converted it for us)
        offset = actualQstart
    else:
        offset = int(qsRec['vqstart'] - qsRec['vstart'])

    offset = max(0, offset)

    vh = record.seq[offset:]

    if len(vh) % 3 != 0:
        vh = vh[:-1 * (len(vh) % 3)]

    if fr4cut and not np.isnan(qsRec['fr4.end']):
        vh = record.seq[offset:int(qsRec['fr4.end'])]

    if primer5seqs:
        primer = str(vh[max(0, end5offset):max(0, end5offset) + maxPrimer5Length])
        qsRec['5endPrimer'], qsRec['5end'], qsRec['5endIndel'] = findBestMatchedPattern(primer, primer5seqs)

    if primer3seqs:
        primer = str(vh[-1 * maxPrimer3Length:])
        qsRec['3endPrimer'], qsRec['3end'], qsRec['3endIndel'] = findBestMatchedPattern(primer, primer3seqs)

    # finish
