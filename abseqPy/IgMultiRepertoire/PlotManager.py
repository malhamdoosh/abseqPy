"""
As of Fri Jun 29 14:42:16 AEST 2018
The only plots that STILL PLOTS IN PYTHON despite pythonPlotOn = False, because abseqR hasn't supported them are:
    1) plotSeqLenDist           - only one that still depends on this is -t seqlen
                                  (abseqR still doesn't have support for this task)
    2) plotVenn                 - used by restriction sites and primer analysis
    4) plotHeatmapFromDF        - restriction sites analysis (detailed and simple)
The plots that OBEY pythonPlotOn() = False is:
    1) all diversity plots (rarefaction, duplication, recapture)
    2) plotDist()
    3) plotSeqLenDistClasses
    4) barLogo                  - generateCumulativeLogo
    5) plotHeatMap              - generateStatsHeatmap
"""


class PlotManager:
    """
    This class acts as a messenger between abseqPy and abseqR.

    It decides whether or not the python backend will be plotting anything (default = no).

    It also has methods that flush the required metadata for abseqR to determine
    what samples are being compared against each other in a file named after the value of _cfg
    """
    _pythonPlotting = False

    def __init__(self):
        PlotManager._pythonPlotting = False

    @staticmethod
    def pythonPlotOn():
        return PlotManager._pythonPlotting

