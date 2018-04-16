import os

from multiprocessing import Queue

from abseq.config import RESULT_FOLDER
from abseq.IgMultiRepertoire.AbSeqWorker import AbSeqWorker, AbSeqWorkerException
from abseq.IgMultiRepertoire.PlotManager import PlotManager
from abseq.IgRepertoire.IgRepertoire import IgRepertoire
from abseq.argsParser import parseYAML, parseArgs


class IgMultiRepertoire:
    def __init__(self, args):
        self.result = Queue()
        self.buffer = []
        self.plotManager = PlotManager(args)
        sampleNames = []
        if args.yaml is not None:
            outdirs = set()
            documents, hasComparisons = parseYAML(args.yaml)
            for yamlArg in documents[:(-1 if hasComparisons else len(documents))]:
                arg = parseArgs(yamlArg)
                arg.outdir = os.path.abspath(arg.outdir) + os.path.sep
                outdirs.add(arg.outdir)
                arg.log = os.path.join(arg.outdir, RESULT_FOLDER, arg.name, arg.name + ".log")
                sampleNames.append(arg.name)
                self.buffer.append(IgRepertoire(**vars(arg)))
            if len(outdirs) == 1:
                outdir = list(outdirs)[0]
            else:
                raise Exception("Multiple output directory in YAML is currently not supported (yet)")
            self.plotManager.processComparisons(documents, sampleNames, hasComparisons, outdir)
        else:
            outdir = args.outdir = os.path.abspath(args.outdir) + os.path.sep
            # <outdir>/result/<sample_name>/<sample_name>.log
            args.log = os.path.join(args.outdir, RESULT_FOLDER, args.name, "{}.log".format(args.name))
            self.buffer.append(IgRepertoire(**vars(args)))
            self.plotManager.processSingleInput(args.name, outdir)
        self.sampleCount = len(self.buffer)

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

        self.result.close()
        self.result.join_thread()

        if noExceptionRaised:
            self.plotManager.plot()

    def rockNRoll(self):

        # initialize workers
        workers = [AbSeqWorker(rep, self.result) for rep in self.buffer]

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
                    raise AbSeqWorkerException(*res)
                self.buffer.append(res)

            for w in workers:
                w.join()
            # done

        except AbSeqWorkerException as e:
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

