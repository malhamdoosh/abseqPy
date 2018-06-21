import os
import ast
import sys

from abseqPy.config import ABSEQROOT

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
    """
    This class acts as a messenger between abseqPy and abseqR.

    It decides whether or not the python backend will be plotting anything (default = no).

    It also has methods that flush the required metadata for abseqR to determine
    what samples are being compared against each other in a file named after the value of _cfg
    """
    # by default, don't plot in python unless rscripting is turned off
    _pythonPlotting = False
    _cfg = "abseq.cfg"

    def __init__(self):
        PlotManager._pythonPlotting = False

    @staticmethod
    def pythonPlotOn():
        return PlotManager._pythonPlotting

    @staticmethod
    def flushSample(name, outdir):
        with open(os.path.join(outdir, PlotManager._cfg), 'w') as fp:
            fp.write(ABSEQROOT + '\n')
            fp.write(name + '\n')

    @staticmethod
    def flushComparisons(pairings, sampleNames, hasComparisons, outdir):
        written = set()
        with open(os.path.join(outdir, PlotManager._cfg), 'w') as fp:
            fp.write(ABSEQROOT + '\n')

            for sample in sampleNames:
                if sample not in written:
                    fp.write(sample + '\n')
                    written.add(sample)

            if hasComparisons:
                if pairings[-1][0] != '--compare':
                    # this will NEVER happen.
                    raise ValueError("Uhhh ... ?")
                for comparison in ast.literal_eval(pairings[-1][1]):
                    userSamples = map(lambda x: x.strip(), comparison.split(","))
                    for s in userSamples:
                        if s not in sampleNames:
                            raise ValueError("Unknown sample name {}, not one of {}".format(s, sampleNames))
                    samples = ','.join(userSamples)
                    if samples not in written:
                        fp.write(samples + '\n')
                        written.add(samples)

