import traceback
import sys

from multiprocessing import Process, current_process


class GeneralWorker(Process):
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

    def __init__(self, jobQueue, resultQueue, jobDescription, *args, **kwargs):
        super(GeneralWorker, self).__init__()
        self.jobQueue = jobQueue
        self.resultQueue = resultQueue
        if jobDescription not in GeneralWorker.ops:
            raise Exception("Unknown job requested {}".format(jobDescription))
        self.jobDescription = jobDescription
        self.args = args
        self.kwargs = kwargs

    def run(self):
        while True:
            job = self.jobQueue.get()
            if job is None:
                break
            try:
                if self.args and not self.kwargs:
                    getattr(job, self.jobDescription)(*self.args)
                elif not self.args and self.kwargs:
                    getattr(job, self.jobDescription)(**self.kwargs)
                else:
                    getattr(job, self.jobDescription)(*self.args, **self.kwargs)
            except Exception as e:
                exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
                fmtMsg = ("Job name: " + str(job.name) + " :: An error occurred while processing " +
                          str(self.jobDescription))
                newE = e.message
                if self.resultQueue is not None:
                    self.resultQueue.put((newE, fmtMsg, traceback.format_exception(exceptionType, exceptionValue,
                                                                                   exceptionTraceback)))
                continue
            # job done
            if self.resultQueue is not None:
                # XXX: cloneAnnot and cloneSeqs are the largest objects in a IgReportoire object,
                # we remove them so that we can pickle them into the queue again.
                # When needed, these files will be loaded automatically later on anyway
                job.cloneAnnot = None
                job.cloneSeqs = None
                self.resultQueue.put(job)
        return


class GeneralWorkerException(Exception):
    def __init__(self, message, errors, tracebackMsg):
        super(GeneralWorkerException, self).__init__(message)
        self.errors = errors
        self.tracebackMsg = ''.join(tracebackMsg)
