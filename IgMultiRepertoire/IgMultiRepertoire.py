from IgRepertoire.IgRepertoire import IgRepertoire
from IgRepertoire.igRepUtils import inferSampleName, detectFileFormat
from multiprocessing import Queue
from copy import deepcopy
from config import DEFAULT_MERGER
from GeneralWorker import GeneralWorker
from math import floor
import os


class IgMultiRepertoire:
    def __init__(self, args):
        self.queue = Queue()
        self.result = Queue()
        self.sampleCount = 0
        self.resource = args.threads
        if os.path.isdir(args.f1):
            allfiles = os.listdir(args.f1)
            clusterFiles = self.__pairFiles(allfiles)
            arg = deepcopy(args)
            arg.merger = args.merger if args.merger is not None else DEFAULT_MERGER
            avgResource = int(floor(self.resource/self.sampleCount))
            avgResource = avgResource if avgResource > 0 else 1
            arg.threads = avgResource
            for sample in clusterFiles:
                self.sampleCount += 1
                if type(sample) == tuple:
                    # paired end sample
                    f1name, f2name = sample
                    retval = inferSampleName(f1name)
                    arg.f1 = f1name
                    arg.f2 = f2name
                    f1Fmt = detectFileFormat(arg.f1)
                    if f1Fmt != detectFileFormat(arg.f2):
                        raise Exception("Detected mismatch in file extensions {} and {}!"
                                        " Both should be either FASTA or FASTQ.".format(arg.f1, arg.f2))
                    arg.fmt = f1Fmt
                else:
                    # single ended
                    f1name = sample
                    retval = inferSampleName(f1name)
                    arg.f1 = f1name
                    arg.f2 = None
                    arg.merger = None
                    arg.fmt = detectFileFormat(f1name)

                arg.outdir += retval[0]
                arg.name = retval[1]
                arg.outdir = (os.path.abspath(args.outdir) + '/').replace("//", "/")
                arg.log = arg.outdir + arg.name + '.log'
                if not os.path.exists(args.outdir):
                    os.makedirs(args.outdir)
                self.queue.put(IgRepertoire(arg))
        else:
            self.sampleCount += 1
            self.queue.put(IgRepertoire(args))

    def analyzeAbundance(self, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.ABUN,  all=all)
        self.__refillQueue()

    def runFastqc(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.FASTQC)
        self.__refillQueue()

    def analyzeDiversity(self, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.DIVER, all=all)
        self.__refillQueue()

    def analyzeProductivity(self, generateReport=True, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.PROD, generateReport=generateReport, all=all)
        self.__refillQueue()

    def annotateClones(self, outDirFilter=None, all=False):
        self.__beginWork(self.queue, self.result, GeneralWorker.ANNOT, outDirFilter=outDirFilter, all=all)
        import sys
        print("Done")
        sys.stdout.flush()
        self.__refillQueue()

    def analyzeRestrictionSites(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.RSA)
        self.__refillQueue()

    def analyzePrimerSpecificity(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.PROD)
        self.__refillQueue()

    def analyze5UTR(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.UTR5)
        self.__refillQueue()

    def analyzeRestrictionSitesSimple(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.RSAS)
        self.__refillQueue()

    def analyzeSecretionSignal(self):
        self.__beginWork(self.queue, self.result, GeneralWorker.SECR)
        self.__refillQueue()

    def __pairFiles(self, files):
        """
        given a list of files, attempt to pair them based on prefix name
        :param files: list of files to be clustered together based on prefixes
        :return: list of files, if element of list is a tuple, then it's detected as paired end,
        or else it will just be a string
        """
        def findPartner(fname):
            lookingFor = "_r2" if '_r1' in fname.lower() else "_r1"
            canonicalName = fname[:fname.lower().rfind("_r1" if lookingFor == '_r2' else '_r2')]
            for f in files:
                if f.lower().startswith(canonicalName) and lookingFor in f.lower():
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

        res = []
        paired = set()
        for f in files:
            if f in paired:     # if this F was already paired with a previous file (paired end), ignore it
                continue
            if "_r1" in f.lower() or "_r2" in f.lower():
                partner = findPartner(f)
                if partner is None:
                    raise Exception("Failed to find opposite read for file {}".format(f))
                res.append(reorderRead(os.path.abspath(f), os.path.abspath(partner)))
                paired.add(partner)
            else:
                # single file
                res.append(os.path.abspath(f))
        return res

    def __beginWork(self, queue, result, jobdesc, *args, **kwargs):
        workers = []
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
            for w in workers:
                import sys
                print("Waiting! Workers: {} worker: {}".format(len(workers), w.name))
                sys.stdout.flush()
                w.join()
                print("DONE! Workers: {}".format(len(workers)))
                sys.stdout.flush()
            # done
        finally:
            for w in workers:
                w.terminate()

    def __refillQueue(self):
        # remove all poison pills from queue
        while not self.queue.empty():
            self.queue.get()

        # sentinel value - None (.empty() and .qsize() not reliable)
        self.__fillPoisonPill(1, self.result)
        itemsAdded = 0  # sanity check
        while True:
            it = self.result.get()
            if it is None:
                break
            self.queue.put(it)
            itemsAdded += 1

        assert itemsAdded == self.sampleCount

    def __fillPoisonPill(self, n, queue, pill=None):
        for _ in range(n):
            queue.put(pill)

