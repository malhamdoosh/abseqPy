import subprocess
import sys
import itertools

from abseq.config import RSCRIPT_SAMPLE_SEPARATOR, ABSEQROOT

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
    _tmpFile = "rscripts_meta.tmp"

    def __init__(self, args):
        self.rscriptArgs = args.rscripts
        # there's no R plotting if user specified that r scripts shouldn't run => python's plotting should run instead
        PlotManager._pythonPlotting = PlotManager.rscriptsOff(self.rscriptArgs)
        self.metadata = []
        self.args = args

    # the following obeys the behaviour of -rs / --rscripts in AbSeq's parser rules.
    # It's trivially easy to write but it improves reading when checking for parser logic
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
        return PlotManager.rscriptsHasArgs(arg) and \
               len(str(arg).split(RSCRIPT_SAMPLE_SEPARATOR)) == 1 and \
               not PlotManager.rscriptsOff(arg)

    @staticmethod
    def rscriptsIsPairedStrings(arg):
        return PlotManager.rscriptsHasArgs(arg) and \
               len(str(arg).split(RSCRIPT_SAMPLE_SEPARATOR)) > 1 and \
               not PlotManager.rscriptsOff(arg)

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
            self._flushMetadata(ABSEQROOT)
            # todo: because main script overridden stdout to AbSeq.log
            # todo: change this to per-sample basis when logging is implemented correctly.
            # i.e. use self.log instead of redirecting to sys.stdout (AbSeq.log)
            retval = subprocess.call(["Rscript",
                                      ABSEQROOT + "/rscripts/masterScript.R"],
                                     stdout=sys.stdout,
                                     stderr=sys.stdout
                                     )
            if retval != 0:
                print("-" * 30)
                print("Error detected in R plotting")
                print("-" * 30)
            # os.remove(PlotManager._tmpFile)        TODO: necessary?

    def processInput(self, allFiles):
        from abseq.IgRepertoire.igRepUtils import inferSampleName
        canonicalNameChangeMap = {}
        for sample in allFiles:
            if type(sample) == tuple:
                # paired end sample
                f1name, _ = sample
                dirName, canonicalName = inferSampleName(f1name, merger=True,
                                                         fastqc=(self.args.task.lower() == 'fastqc'))
                outDirName = self.args.outdir + dirName
                self.metadata.append((outDirName, canonicalName))
            else:
                dirName, canonicalName = inferSampleName(sample, merger=False,
                                                         fastqc=(self.args.task.lower() == 'fastqc'))
                outDirName = self.args.outdir + dirName
                self.metadata.append((outDirName, canonicalName))
            canonicalNameChangeMap[canonicalName] = canonicalName

        # 1st. Filter out all samples that are not required. (by consulting -rs <args>)
        # 2nd. Post filtering, remap all abseq canonical names to user provided names
        if self.rscriptArgs and self.rscriptArgs != 'off':
            # take away samples that are not requested by user if -rs was specified
            requestedSamples = self._getRscriptSamples()
            self.metadata = filter(lambda x: x[1] in requestedSamples, self.metadata)
            canonicalNameChangeMap = {k: v for k, v in canonicalNameChangeMap.items() if k in requestedSamples}

            # remap abseq's infered names to user provided ones - flatten rscriptArgs
            for userProvidedNames in list(itertools.chain(*self.rscriptArgs)):
                res = self._findBestMatch(userProvidedNames)[1]
                # update new name
                newName = self._findBestMatch(userProvidedNames, useProvidedName=True)[1]
                # if user provided confusing names, eg: PCR1_L001 and PCR1 both refer to the same sample,
                # we raise an exception immediately!
                if res in canonicalNameChangeMap and canonicalNameChangeMap[res] != res and \
                        canonicalNameChangeMap[res] != newName:
                    raise Exception("Misleading -rs argument, {} and {}"
                                    .format(self._findBestMatch(userProvidedNames, useProvidedName=True)[1],
                                            canonicalNameChangeMap[res]))
                canonicalNameChangeMap[res] = newName

            self.metadata = map(lambda x: (x[0], canonicalNameChangeMap[x[1]]), self.metadata)

        return canonicalNameChangeMap

    def addMetadata(self, data):
        self.metadata.append(data)

    def _getRscriptSamples(self):
        # whether or not there was python plotting, see if user explicitly chose samples
        requestedSamples = set()
        if self.rscriptArgs and self.rscriptArgs != 'off':
            for pairings in self.rscriptArgs:
                for sampleName in pairings:
                    res = self._findBestMatch(sampleName)
                    if res:
                        requestedSamples.add(res[1])
        return requestedSamples

    def _flushMetadata(self, abSeqRootDir):
        """
        AbSeq's way of communicating with external Rscript. AbSeq's python component supplies R with
        all the necessary information to find (CSVs) and generate plots.

        The first line is the absolute path to abseq's root directory.

        Output format is of:
        /path/to/sampleDirectory1, /path/to/sampleDirectory2, ... /path/to/sampleDirectoryN ? sampleName1, ..., sampleNameN

        Where commas separate pairings of samples and question mark denotes that everything on RHS is the
        canonical sample name used throughout AbSeq for each given sample (in the order they appeared on LHS)

        Non-paired samples (even paired samples requires R to plot sample individually) will have the format of:

        /path/to/sampleDirectory/ ? sampleName
        (it's the same as multi-sample except that there are no commas)
        :return: None. Side effect : writes a temporary file
        """
        # only write metadata if we have to plot in R
        if not self._pythonPlotting:
            writeBuffer = []
            # refine the entries provided by user in self.rscripts' argument to the canonical name
            # as defined in AbSeq and the directory this sample lives in
            if self.rscriptArgs:
                for pairings in self.rscriptArgs:
                    # this if statement will be false when someone purposely provided a
                    # SINGLE sample as "comparison/pairing" E.G. -rs "PCR1 |" or -rs "PCR1 | ; PCR2" etc..
                    if not (len(pairings) == 1 or '' in pairings):
                        writeBuffer.append([self._findBestMatch(sampleName, useProvidedName=True)
                                            for sampleName in pairings])

                # at this point, writeBuffer is a list of list as such:
                # writeBuffer = [
                #           [ (self.outdir/PCR1_BZ123_ACGGCT_GCGTA_L001, PCR1_L001),
                #               (self.outdir/PCR2_BZC1_ACGGTA_GAGA_L001, PCR2_L001), .. ]
                #           [ (self.outdir/PCR4_BZ123_ACGG_ACGG_L001, PCR4_L001),
                #               (self.outdir/PCR5_....., PCR5_L001), ... ]
                #           [ ... ]
                # ] - each inner list is a set of pairing

            # the individual samples has to be plotted too!
            writeBuffer += map(lambda x: [x], self.metadata)

            with open(PlotManager._tmpFile, "w") as fp:
                # tell R where AbSeq lives
                fp.write(abSeqRootDir + "\n")
                for differentPairings in writeBuffer:
                    # write all directories for a given pairing, then the canonical name, separated by a '?' token
                    fp.write(','.join(map(lambda x: x[0].lstrip("/") if x else '', differentPairings)) + "?")
                    fp.write(','.join(map(lambda x: x[1] if x else '', differentPairings)) + "\n")
                    # final result, rscripts_meta.tmp looks like:
                    # (pairings come first)
                    # self.outdir/PCR1_BZ123_ACGGCT_GCGTA_L001,self.outdir/PCR2_BZC1_ACGGTA_GAGA_L001,...?PCR1_L001,...
                    # self.outdir/PCR4_BZ123_ACGG_ACGG_L001,self.outdir/PCR5_...,..?PCR4_L001,PCR5_L001,...
                    # .
                    # .
                    # .
                    # (then single samples)
                    # self.outdir/PCR1_BZ123_..._L001,PCR1_L001
                    # self.outdir/PCR2_BZ123_..._L001,PCR2_L001

    def _findBestMatch(self, sampleName, useProvidedName=False):
        v = float('-inf')
        bestMatch = None
        if sampleName:
            for sampleDir, sampleCName in self.metadata:
                matchScore = max(_nameMatch(sampleDir, sampleName), _nameMatch(sampleCName, sampleName))
                if matchScore > v:
                    bestMatch = (sampleDir, sampleName if useProvidedName else sampleCName)
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
            maxval = max(maxval, matrix[i][j])
    return maxval



