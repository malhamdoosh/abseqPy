from PlotManager import PlotManager
from IgRepertoire.IgRepertoire import IgRepertoire
from IgRepertoire.igRepUtils import inferSampleName, detectFileFormat
from multiprocessing import Queue
from copy import deepcopy
from config import DEFAULT_MERGER
from GeneralWorker import GeneralWorker, GeneralWorkerException
from math import floor
import os


class IgMultiRepertoire:
    def __init__(self, args):
        self.queue = Queue()
        self.result = Queue()
        self.sampleCount = 0
        self.resource = args.threads
        self.plotManager = PlotManager(args)
        if os.path.isdir(args.f1):
            clusterFiles = self.__pairFiles(args.f1, args)
            self.sampleCount = len(clusterFiles)
            avgResource = int(floor(self.resource / self.sampleCount))
            avgResource = avgResource if avgResource > 0 else 1
            for sample in clusterFiles:
                modifiedArgs = deepcopy(args)
                modifiedArgs.merger = args.merger if args.merger is not None else DEFAULT_MERGER
                modifiedArgs.threads = avgResource
                if type(sample) == tuple:
                    # paired end sample
                    f1name, f2name = sample
                    retval = inferSampleName(f1name, merger=True, fastqc=(args.task.lower() == 'fastqc'))
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
                    retval = inferSampleName(f1name, merger=False, fastqc=(args.task.lower() == 'fastqc'))
                    modifiedArgs.f1 = f1name
                    modifiedArgs.f2 = None
                    modifiedArgs.merger = None
                    modifiedArgs.fmt = detectFileFormat(f1name)
                modifiedArgs.outdir += retval[0]
                modifiedArgs.name = retval[1]
                modifiedArgs.outdir = (os.path.abspath(modifiedArgs.outdir) + '/').replace("//", "/")
                modifiedArgs.log = modifiedArgs.outdir + modifiedArgs.name + '.log'
                self.plotManager.addMetadata(retval)
                if not os.path.exists(modifiedArgs.outdir):
                    os.makedirs(modifiedArgs.outdir)
                self.queue.put(IgRepertoire(modifiedArgs))
        else:
            self.plotManager.addMetadata(
                (inferSampleName(args.f1, args.merger, args.task.lower() == 'fastqc')[0], args.name))
            self.sampleCount += 1
            self.queue.put(IgRepertoire(args))

    def analyzeAbundance(self, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.ABUN, all=all)

    def runFastqc(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.FASTQC)

    def analyzeDiversity(self, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.DIVER, all=all)

    def analyzeProductivity(self, generateReport=True, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.PROD, generateReport=generateReport, all=all)

    def annotateClones(self, outDirFilter=None, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.ANNOT, outDirFilter=outDirFilter, all=all)

    def analyzeRestrictionSites(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.RSA)

    def analyzePrimerSpecificity(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.PRIM)

    def analyze5UTR(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.UTR5)

    def analyzeRestrictionSitesSimple(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.RSAS)

    def analyzeSecretionSignal(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.SECR)

    def analyzeSeqLen(self, klass=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.SEQLEN, klass=klass)

    def finish(self):
        """
        Queue might still be buffered, finish off cleanup here.
        Then, delegate to plot manager to decide if there's further plotting required.
        :return: None
        """
        # pop all items
        while not self.queue.empty():
            self.queue.get()
        while not self.result.empty():
            self.result.get()
        self.queue.close()
        self.queue.join_thread()
        self.result.close()
        self.result.join_thread()
        self.plotManager.plot()

    def __pairFiles(self, folder, args):
        """
        given a list of files, attempt to pair them based on prefix name
        :param folder: folder in which these files are found in
        :return: list of files, if element of list is a tuple, then it's detected as paired end,
        or else it will just be a string
        """
        files = os.listdir(os.path.abspath(folder))

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
                res.append(reorderRead(os.path.abspath(folder + "/" + f), os.path.abspath(folder + "/" + partner)))
                paired.add(partner)
            else:
                # single file
                res.append(os.path.abspath(folder + '/' + f))

        return distinct(res)

    def __beginWork(self, queue, result, jobdesc, refill=True, *args, **kwargs):
        workers = []
        collectedResult = []
        # initialize workers
        for _ in range(min(self.sampleCount, self.resource)):
            workers.append(GeneralWorker(queue, result, jobdesc, *args, **kwargs))

        # since macOSX doesn't support .qsize(). also, .empty() and .qsize() are
        # unreliable ==> we use poison pills
        self.__fillPoisonPill(self.resource + 10, self.queue)

        try:
            # start workers
            for w in workers:
                w.start()

            # wait for all workers to complete
            for i in xrange(self.sampleCount):
                res = self.result.get()
                if type(res) == tuple:
                    # XXX: encountered an exception! - here, decide to raise it immediately.
                    # all accompanying processes will halt immediately due to this raise.
                    raise GeneralWorkerException(*res)
                if refill:
                    collectedResult.append(res)

            for w in workers:
                w.join()
            # done

            # refill if needed
            if refill:
                # empty out original queue from Nones
                while not self.queue.empty():
                    res = self.queue.get()
                    assert res is None

                itemsAdded = 0  # sanity check
                for i in collectedResult:
                    self.queue.put(i)
                    itemsAdded += 1

                assert itemsAdded == self.sampleCount
        except GeneralWorkerException as e:
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

    def __fillPoisonPill(self, n, queue, pill=None):
        for _ in range(n):
            queue.put(pill)
