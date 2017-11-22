import os


"""
XXX: IMPORTANT NOTE
As of Nov 22 2017: 
The only plots that STILL PLOTS IN PYTHON despite pythonPlotOn = False is:
    1) plotSeqLenDist           -- will be moved to R #TODO
    2) plotSeqLenDistClasses    -- will be moved to R #TODO
    3) plotVenn                 - mostly used by restriction sites
    4) plotHeatMap              - generateStatsHeatmap
    5) plotHeatmapFromDF
    6) barLogo                  - generateCumulativeLogo
The plots that OBEY pythonPlotOn() = False is:
    1) all diversity plots (rarefaction, duplication, recapture)
    2) plotDist()
"""


class PlotManager:
    # by default, don't plot in python unless rscripting is turned off
    _pythonPlotting = False

    def __init__(self, rscriptArgs):
        self.rscriptArgs = rscriptArgs
        # there's no R plotting if user specified that r scripts shouldn't run => python's plotting should run instead
        PlotManager._pythonPlotting = PlotManager.rscriptsOff(self.rscriptArgs)
        self.metadata = []

    # the following obeys the behaviour of -rs / --rscripts in AbSeq's parser rules.
    # It's trivially easy to write but it improves reading when checking for parser logic above.
    @staticmethod
    def rscriptsHasNoArgs(arg):
        return arg is None

    @staticmethod
    def rscriptsHasArgs(arg):
        return arg is not None

    @staticmethod
    def rscriptsOff(arg):
        return PlotManager.rscriptsHasArgs(arg) and str(arg).lower() == 'off'

    @staticmethod
    def rscriptsIsConf(arg):
        return PlotManager.rscriptsHasArgs(arg) and len(str(arg).split(",")) == 1 and os.path.exists(arg)

    @staticmethod
    def rscriptsIsPairedStrings(arg):
        return PlotManager.rscriptsHasArgs(arg) and len(str(arg).split(",")) > 1

    @staticmethod
    def pythonPlotOn():
        """
        multiple threads / process may read this value. It's a RO value - it'll never be changed after RScriptsManager
        is instantiated. There's no need to guard this against a race condition.
        :return: True if python should plot graphs, false otherwise.
        """
        return PlotManager._pythonPlotting

    def addMetadata(self, sample):
        """
        stores metadata about a given sample. More specifically, stores the output directory
        and canonical name of a sample so that the rscripts that will plot multiple-repertoire
        comparison will know which (sample) and where to find them
        :param sample: sample is a tuple type of (sample's directory, sample's canonical name)
        :return: None
        """
        # if not self.noPlot:
        self.metadata.append(sample)

    def flushMetadata(self):
        # only write metadata if we have to plot in R
        if not self._pythonPlotting and type(self.rscriptArgs) == list:
            writeBuffer = []
            # refine the entries provided by user in self.rscripts' argument to the canonical name
            # as defined in AbSeq and the directory this sample lives in
            for pairings in self.rscriptArgs:
                writeBuffer.append([self.__findBestMatch(sampleName) for sampleName in pairings])
            # at this point, writeBuffer is a list of list as such:
            # writeBuffer = [
            #           [ (PCR1_BZ123_ACGGCT_GCGTA_L001/, PCR1_L001), (PCR2_BZC1_ACGGTA_GAGA_L001/, PCR2_L001), .. ]
            #           [ (PCR4_BZ123_ACGG_ACGG_L001/, PCR4_L001), (PCR5_....., PCR5_L001), ... ]
            #           [ ... ]
            # ] - each inner list is a set of pairing
            with open("rscripts_meta.tmp", "w") as fp:
                for pairing in writeBuffer:
                    # write all directories for a given pairing, then the canonical name, separated by a '?' token
                    fp.write(','.join(map(lambda x: x[0], pairing)) + "?")
                    fp.write(','.join(map(lambda x: x[1], pairing)) + "\n")

                    # final result, rscripts_meta.tmp looks like:
                    # PCR1_BZ123_ACGGCT_GCGTA_L001/, PCR2_BZC1_ACGGTA_GAGA_L001/, ... ? PCR1_L001, PCR2_L001, ...
                    # PCR4_BZ123_ACGG_ACGG_L001/, PCR5_.../, .. ? PCR4_L001, PCR5_L001, ...
                    # ...
    def __findBestMatch(self, sampleName):
        v = float('-inf')
        bestMatch = None
        for sampleDir, sampleCName in self.metadata:
            matchScore = max(_nameMatch(sampleDir, sampleName), _nameMatch(sampleCName, sampleName))
            if matchScore > v:
                bestMatch = (sampleDir, sampleCName)
                v = matchScore
        return bestMatch


def _nameMatch(string1, string2, deletionPenalty=-2, insertionPenalty=-2, matchScore=5, mismatchScore=-3):
    """
    returns local edit distance between 2 strings
    :param string1: a string type
    :param string2: a string type
    :param deletionPenalty: penalty score for a deletion
    :param insertionPenalty: penalty score for an insertion
    :param matchScore: reward score for each match
    :param mismatchScore: penalty score for a substitution
    :return: Local edit distance between string 1 and string 2
    """
    string1 = string1.strip()
    string2 = string2.strip()
    matrix = [[0] * (len(string1) + 1)] * (len(string2) + 1)
    maxval = 0
    for i in xrange(1, len(string2) + 1):
        for j in xrange(1, len(string1) + 1):
            matrix[i][j] = max(matrix[i][j - 1] + deletionPenalty, matrix[i - 1][j] + insertionPenalty,
                               matrix[i - 1][j - 1] + (
                                   matchScore if string1[j - 1] == string2[i - 1] else mismatchScore))
            if matrix[i][j] > maxval:
                maxval = matrix[i][j]
    return maxval
