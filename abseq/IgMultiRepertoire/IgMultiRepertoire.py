import os

from multiprocessing import Queue
from copy import deepcopy
from math import floor

from abseq.config import DEFAULT_MERGER
from abseq.IgMultiRepertoire.GeneralWorker import GeneralWorker, GeneralWorkerException
from abseq.IgMultiRepertoire.PlotManager import PlotManager
from abseq.IgRepertoire.IgRepertoire import IgRepertoire
from abseq.IgRepertoire.igRepUtils import inferSampleName, detectFileFormat


class IgMultiRepertoire:
    def __init__(self, args):
        self.queue = Queue()
        self.result = Queue()
        self.buffer = []
        self.sampleCount = 0
        self.resource = args.threads
        self.plotManager = PlotManager(args)

        if os.path.isdir(args.f1):
            clusterFiles = self._pairFiles(args.f1, args)
            canonicalNameChangeMap = self.plotManager.processInput(clusterFiles)

            # get requested samples from -rs (if specified / if any) only
            self.sampleCount = len(canonicalNameChangeMap)

            # calculate resource (averaged over all samples)
            avgResource = int(floor(self.resource / self.sampleCount))
            avgResource = avgResource if avgResource > 0 else 1
            for sample in clusterFiles:
                modifiedArgs = deepcopy(args)
                modifiedArgs.merger = args.merger if args.merger is not None else DEFAULT_MERGER
                modifiedArgs.threads = avgResource
                if type(sample) == tuple:
                    # paired end sample
                    f1name, f2name = sample
                    inferredDir, inferredName = inferSampleName(f1name, merger=True,
                                                                fastqc=(args.task.lower() == 'fastqc'))
                    modifiedArgs.f1 = f1name
                    modifiedArgs.f2 = f2name
                    f1Fmt = detectFileFormat(modifiedArgs.f1)
                    if f1Fmt != detectFileFormat(modifiedArgs.f2):
                        raise Exception("Detected mismatch in file extensions {} and {}!"
                                        " Both should be either FASTA or FASTQ.".format(modifiedArgs.f1,
                                                                                        modifiedArgs.f2))
                    modifiedArgs.fmt = f1Fmt
                else:
                    # single ended
                    f1name = sample
                    inferredDir, inferredName = inferSampleName(f1name, merger=False,
                                                                fastqc=(args.task.lower() == 'fastqc'))
                    modifiedArgs.f1 = f1name
                    modifiedArgs.f2 = None
                    modifiedArgs.merger = None
                    modifiedArgs.fmt = detectFileFormat(f1name)

                # if abseq's inferred sample name is in the map ==> sample was specified in -rs
                # we also need to remap the name
                if inferredName in canonicalNameChangeMap:
                    modifiedArgs.outdir += inferredDir
                    modifiedArgs.name = canonicalNameChangeMap[inferredName]
                    modifiedArgs.outdir = os.path.abspath(modifiedArgs.outdir) + os.path.sep
                    modifiedArgs.log = modifiedArgs.outdir + modifiedArgs.name + '.log'
                    self.buffer.append(IgRepertoire(**vars(modifiedArgs)))
        else:
            self.plotManager.addMetadata((args.outdir, args.name))
            args.outdir = os.path.abspath(args.outdir) + os.path.sep
            args.log = args.outdir + args.name + ".log"
            self.sampleCount += 1
            self.buffer.append(IgRepertoire(**vars(args)))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Extremely important that finish is called (for 3 reasons):
            1. Make sure that primer specificity was conducted if either one of primer files were provided (regardless of
               the -t <task> option) - unless, of course, if -t was already 'primer',
               then no additional analysis is required.
            2. Queue might still be buffered, finish off cleanup here.
            3. Then, delegate to plot manager to decide if there's further plotting required.
        :return: None
        """
        noExceptionRaised = exc_tb is None and exc_val is None and exc_tb is None

        if noExceptionRaised:
            # make sure that if user specified either one of primer end file, we unconditionally run primer analysis
            # (duh)
            if (self.buffer[0].end3 or self.buffer[0].end5) and self.buffer[0].task != 'primer':
                print("Primer file detected, conducting primer specificity analysis ... ")
                self.analyzePrimerSpecificity()

        self.queue.close()
        self.queue.join_thread()
        self.result.close()
        self.result.join_thread()

        if noExceptionRaised:
            self.plotManager.plot()

    def analyzeAbundance(self):
        self._beginWork(GeneralWorker.ABUN)

    def runFastqc(self):
        self._beginWork(GeneralWorker.FASTQC)

    def analyzeDiversity(self):
        self._beginWork(GeneralWorker.DIVER)

    def analyzeProductivity(self):
        self._beginWork(GeneralWorker.PROD)

    def annotateClones(self, outDirFilter=None):
        self._beginWork(GeneralWorker.ANNOT, outDirFilter=outDirFilter)

    def analyzeRestrictionSites(self):
        self._beginWork(GeneralWorker.RSA)

    def analyzePrimerSpecificity(self):
        self._beginWork(GeneralWorker.PRIM)

    def analyze5UTR(self):
        self._beginWork(GeneralWorker.UTR5)

    def analyzeRestrictionSitesSimple(self):
        self._beginWork(GeneralWorker.RSAS)

    def analyzeSecretionSignal(self):
        self._beginWork(GeneralWorker.SECR)

    def analyzeSeqLen(self, klass=False):
        self._beginWork(GeneralWorker.SEQLEN, klass=klass)

    def finish(self):
        self.__exit__(None, None, None)

    def _pairFiles(self, folder, args):
        """
        given a list of files, attempt to pair them based on prefix name
        :param folder: folder in which these files are found in
        :return: list of files, if element of list is a tuple, then it's detected as paired end,
        or else it will just be a string
        """
        files = [f for f in os.listdir(os.path.abspath(folder)) if not f.startswith(".")]

        def findPartner(fname):
            lookingFor = "_r2" if '_r1' in fname.lower() else "_r1"
            canonicalName = fname[:fname.lower().rfind("_r1" if lookingFor == '_r2' else '_r2')]
            for f in files:
                if f.lower().startswith(canonicalName.lower()) and lookingFor in f.lower():
                    return f
            return None

        def reorderRead(a, b):
            """
            return reads as (_R1, _R2)
            :param a: arbitrary file (R1/R2)
            :param b: arbitrary file (R1/R2)
            :return: Ordered file (R1, R2)
            """
            if '_r1' in a.lower():
                return a, b
            return b, a

        def distinct(files):
            seen = set()
            reduced = []
            for f in files:
                if type(f) == tuple:
                    _, sampleName = inferSampleName(f[0], args.merger, args.task.lower() == 'fastqc')
                else:
                    _, sampleName = inferSampleName(f, args.merger, args.task.lower() == 'fastqc')
                if sampleName not in seen:
                    seen.add(sampleName)
                    reduced.append(f)
            return reduced

        def acceptedFormat(filename):
            return detectFileFormat(filename, noRaise=True) is not None

        res = []
        paired = set()
        for f in files:
            if f in paired or not acceptedFormat(f):
                continue
            if "_r1" in f.lower() or "_r2" in f.lower():
                partner = findPartner(f)
                if partner is None:
                    raise Exception("Failed to find opposite read for file {}".format(f))
                res.append(reorderRead(os.path.abspath(os.path.join(folder,  f)),
                                       os.path.abspath(os.path.join(folder, partner))))
                paired.add(partner)
            else:
                # single file
                res.append(os.path.abspath(os.path.join(folder, f)))

        return distinct(res)

    def _beginWork(self, jobdesc, *args, **kwargs):

        # fill self.queue with data from self.buffer
        assert self.queue.empty()
        assert len(self.buffer) == self.sampleCount
        for _ in range(len(self.buffer)):
            self.queue.put(self.buffer.pop())
        assert len(self.buffer) == 0

        # initialize workers
        workers = [GeneralWorker(self.queue, self.result, jobdesc, *args, **kwargs) for _ in
                   range(min(self.sampleCount, self.resource))]

        # since macOSX doesn't support .qsize(). also, .empty() and .qsize() are
        # unreliable ==> we use poison pills
        self._fillPoisonPill(self.resource + 10, self.queue)

        try:
            # start workers
            for w in workers:
                w.start()

            # wait for all workers to complete
            for i in range(self.sampleCount):
                res = self.result.get()
                if type(res) == tuple:
                    # XXX: encountered an exception! - here, decide to raise it immediately.
                    # all accompanying processes will halt immediately due to this raise.
                    raise GeneralWorkerException(*res)
                self.buffer.append(res)

            for w in workers:
                w.join()
            # done

            # empty out original queue from Nones
            while not self.queue.empty():
                res = self.queue.get()
                assert res is None

        except GeneralWorkerException as e:
            print("\n\n{}".format(e.errors))
            print("\n\nSomething went horribly wrong while trying to run AbSeq!")
            print("GeneralWorker stacktrace:")
            print("-" * 120)
            print(e.tracebackMsg)
            print("-" * 120)
            # re-raise exception
            raise e
        # catch-all exception
        except Exception as e:
            raise e
        finally:
            for w in workers:
                w.terminate()

    @staticmethod
    def _fillPoisonPill(n, queue, pill=None):
        for _ in range(n):
            queue.put(pill)
