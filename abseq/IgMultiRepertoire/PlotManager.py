import os
import ast
import sys

from abseq.config import ABSEQROOT

"""
XXX: IMPORTANT NOTE
As of Nov 22 2017: 
The only plots that STILL PLOTS IN PYTHON despite pythonPlotOn = False is:
    1) plotSeqLenDist
    2) plotVenn                 - mostly used by restriction sites
    3) plotHeatMap              - generateStatsHeatmap
    4) plotHeatmapFromDF
The plots that OBEY pythonPlotOn() = False is:
    1) all diversity plots (rarefaction, duplication, recapture)
    2) plotDist()
    3) plotSeqLenDistClasses
    4) barLogo                  - generateCumulativeLogo
"""


class PlotManager:
    # by default, don't plot in python unless rscripting is turned off
    _pythonPlotting = False
    _cfg = "abseq.cfg"

    def __init__(self, args):
        PlotManager._pythonPlotting = False
        self.metadata = []
        self.end5File = args.primer5end
        self.end3File = args.primer3end
        self.upstream = args.upstream
        self.outputDir = args.outdir
        self.args = args

    @staticmethod
    def pythonPlotOn():
        """
        multiple threads / process may read this value. It's a RO value - it'll never be changed after RScriptsManager
        is instantiated. There's no need to guard this against a race condition.
        :return: True if python should plot graphs, false otherwise.
        """
        return PlotManager._pythonPlotting

    def plot(self):
        """
        plots csv files. If python plotting was on, then this function will have no effect
        :return: None. Side effects: plots in R if specified as such
        """
        if not PlotManager._pythonPlotting:
            # todo: because main script overridden stdout to AbSeq.log
            # todo: change this to per-sample basis when logging is implemented correctly.
            # i.e. use self.log instead of redirecting to sys.stdout (AbSeq.log)
            DUMMY_VALUE = "None"
            # primer args
            arg1 = str(self.end5File)
            arg2 = str(self.end3File)

            # upstream args
            arg3 = "{:.0f}".format(self.upstream[0]) if self.upstream else DUMMY_VALUE
            arg4 = "{:.0f}".format(self.upstream[1]) if self.upstream else DUMMY_VALUE
            sys.stdout.flush()
            # retval = subprocess.call(["Rscript",
            #                           ABSEQROOT + "/rscripts/masterScript.R",
            #                           arg1,
            #                           arg2,
            #                           arg3,
            #                           arg4],
            #                          stdout=sys.stdout,
            #                          stderr=sys.stdout
            #                          )
            # if retval != 0:
            #     print("-" * 30)
            #     print("Error detected in R plotting")
            #     print("-" * 30)
            # os.remove(PlotManager._tmpFile)        TODO: necessary?

    def processSingleInput(self, name, outdir):
        with open(os.path.join(outdir, PlotManager._cfg), 'w') as fp:
            fp.write(ABSEQROOT + '\n')
            fp.write(name + '\n')

    def processComparisons(self, pairings, sampleNames, hasComparisons, outdir):
        with open(os.path.join(outdir, PlotManager._cfg), 'w') as fp:
            fp.write(ABSEQROOT + '\n')

            if hasComparisons:
                if pairings[-1][0] != '--compare':
                    # this will NEVER happen.
                    raise ValueError("Uhhh ... ?")
                for comparison in ast.literal_eval(pairings[-1][1]):
                    userSamples = map(lambda x: x.strip(), comparison.split(","))
                    for s in userSamples:
                        if s not in sampleNames:
                            raise ValueError("Unknown sample name {}, not one of {}".format(s, sampleNames))
                    fp.write(','.join(userSamples) + '\n')

            for sample in sampleNames:
                fp.write(sample + '\n')

