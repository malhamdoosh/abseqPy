from multiprocessing import Process, current_process


class GeneralWorker(Process):
    ops = FASTQC, ANNOT, ABUN, PROD, DIVER, SECR, UTR5, RSAS, RSA, PRIM, SEQLEN = 'runFastqc', 'annotateClones', \
                                                                          'analyzeAbundance',\
                                                                          'analyzeProductivity',\
                                                                          'analyzeDiversity', \
                                                                          'analyzeSecretionSignal',\
                                                                          'analyze5UTR', \
                                                                          'analyzeRestrictionSitesSimple', \
                                                                          'analyzePrimerSpecificity', \
                                                                          'analyzeRestrictionSites', \
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
                print("Got job: {} args: {} kwargs: {}".format(self.jobDescription, self.args, self.kwargs))
                import sys
                sys.stdout.flush()
                if self.args and not self.kwargs:
                    getattr(job, self.jobDescription)(*self.args)
                elif not self.args and self.kwargs:
                    getattr(job, self.jobDescription)(**self.kwargs)
                else:
                    getattr(job, self.jobDescription)(*self.args, **self.kwargs)
            except Exception as e:
                print(str(job.name) + ":: An error occurred while processing " + str(self.jobDescription))
                print("Error description: {}".format(str(e)))
                if self.resultQueue is not None:
                    self.resultQueue.put(None)
                continue
            # job done
            if self.resultQueue is not None:
                self.resultQueue.put(job)
        return
