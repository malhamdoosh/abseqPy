from __future__ import print_function
import os
import sys

from multiprocessing import Queue

from abseqPy.config import AUX_FOLDER
from abseqPy.IgMultiRepertoire.AbSeqWorker import AbSeqWorker, AbSeqWorkerException, ResourcePool
from abseqPy.IgRepertoire.IgRepertoire import IgRepertoire
from abseqPy.argsParser import parseYAML, parseArgs


class IgMultiRepertoire:
    def __init__(self, args):
        self.result = Queue()
        self.buffer = []
        sampleNames = []
        if args.yaml is not None:
            outdirs = set()
            for yamlArg in parseYAML(args.yaml):
                arg = parseArgs(yamlArg)
                arg.outdir = os.path.abspath(arg.outdir) + os.path.sep
                outdirs.add(arg.outdir)
                arg.log = os.path.join(arg.outdir, AUX_FOLDER, arg.name, arg.name + ".log")
                sampleNames.append(arg.name)
                self.buffer.append(IgRepertoire(**vars(arg)))
        else:
            # <outdir>/result/<sample_name>/<sample_name>.log
            args.log = os.path.join(args.outdir, AUX_FOLDER, args.name, "{}.log".format(args.name))
            self.buffer.append(IgRepertoire(**vars(args)))
        self.sampleCount = len(self.buffer)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # noExceptionRaised = exc_tb is None and exc_val is None and exc_tb is None
        self.result.close()
        self.result.join_thread()

    def start(self):

        # resource pool - initially all consumed by repertoire objects
        resourcePool = ResourcePool(0)

        # initialize workers
        workers = [AbSeqWorker(rep, self.result, resourcePool) for rep in self.buffer]

        try:
            # start workers
            for w in workers:
                w.start()

            # wait for all workers to complete
            for i in range(self.sampleCount):
                res = self.result.get()
                if isinstance(res, tuple):
                    # XXX: encountered an exception! - here, decide to raise it immediately.
                    # all accompanying processes will halt immediately due to this raise.
                    raise AbSeqWorkerException(*res)
                self.buffer.append(res)

            for w in workers:
                w.join()
            # done

        except AbSeqWorkerException as e:
            print("\n\n{}".format(e.errors), file=sys.stderr)
            print("\n\nSomething went horribly wrong while trying to run AbSeq!", file=sys.stderr)
            print("GeneralWorker stacktrace:", file=sys.stderr)
            print("-" * 120, file=sys.stderr)
            print(e.tracebackMsg, file=sys.stderr)
            print("-" * 120, file=sys.stderr)
            # re-raise exception
            raise e
        # catch-all exception
        except Exception as e:
            raise e
        finally:
            for w in workers:
                w.terminate()

