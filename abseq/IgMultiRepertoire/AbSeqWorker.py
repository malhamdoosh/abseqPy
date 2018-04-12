import traceback
import sys

from multiprocessing import Process, current_process


class AbSeqWorker(Process):
    ops = FASTQC, ANNOT, ABUN, PROD, DIVER, SECR, UTR5, RSAS, RSA, PRIM, SEQLEN = 'runFastqc', 'annotateClones', \
                                                                          'analyzeAbundance',\
                                                                          'analyzeProductivity',\
                                                                          'analyzeDiversity', \
                                                                          'analyzeSecretionSignal',\
                                                                          'analyze5UTR', \
                                                                          'analyzeRestrictionSitesSimple', \
                                                                          'analyzeRestrictionSites', \
                                                                          'analyzePrimerSpecificity', \
                                                                          'analyzeSeqLen'

    def __init__(self, repertoire, resultQueue):
        super(AbSeqWorker, self).__init__()
        self.repertoire = repertoire
        self.resultQueue = resultQueue

    def run(self):
        while True:

            task, args, kwargs = self.repertoire._nextTask()
            if task is None:
                # all jobs done for this repertoire
                self.repertoire._minimize()
                self.resultQueue.put(self.repertoire)
                break

            try:
                getattr(self.repertoire, task)(*args, **kwargs)
            except Exception as e:
                exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
                fmtMsg = ("Job name: " + str(self.repertoire.name) + " :: An error occurred while processing " +
                          str(task))
                newE = e.message
                if self.resultQueue is not None:
                    self.resultQueue.put((newE, fmtMsg, traceback.format_exception(exceptionType, exceptionValue,
                                                                                   exceptionTraceback)))
                continue
        return


class AbSeqWorkerException(Exception):
    def __init__(self, message, errors, tracebackMsg):
        super(AbSeqWorkerException, self).__init__(message)
        self.errors = errors
        self.tracebackMsg = ''.join(tracebackMsg)
