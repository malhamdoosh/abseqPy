import subprocess
import sys
import itertools
import numpy as np
import os

from collections import defaultdict

from abseq.config import RSCRIPT_SAMPLE_SEPARATOR, ABSEQROOT, RSCRIPT_PAIRING_SEPARATOR, RESULT_FOLDER

"""
XXX: IMPORTANT NOTE
As of Nov 22 2017: 
The only plots that STILL PLOTS IN PYTHON despite pythonPlotOn = False is:
    1) plotSeqLenDist           -- will be moved to R #TODO
    3) plotVenn                 - mostly used by restriction sites
    4) plotHeatMap              - generateStatsHeatmap
    5) plotHeatmapFromDF
    6) barLogo                  - generateCumulativeLogo
The plots that OBEY pythonPlotOn() = False is:
    1) all diversity plots (rarefaction, duplication, recapture)
    2) plotDist()
    3) plotSeqLenDistClasses
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
        self.end5File = args.primer5end
        self.end3File = args.primer3end
        self.upstream = args.upstream
        self.args = args
        self._nameToFileMap = {}

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
               len(str(arg).split(RSCRIPT_PAIRING_SEPARATOR)) == 1 and \
               len(str(arg).split(RSCRIPT_SAMPLE_SEPARATOR)) == 1 and \
               not PlotManager.rscriptsOff(arg)

    @staticmethod
    def rscriptsIsPairedStrings(arg):
        sarg = str(arg)
        return PlotManager.rscriptsHasArgs(arg) \
               and (len(sarg.split(RSCRIPT_SAMPLE_SEPARATOR)) > 1 or len(sarg.split(RSCRIPT_PAIRING_SEPARATOR)) > 1) \
               and not PlotManager.rscriptsOff(arg)

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
            DUMMY_VALUE = "None"
            # primer args
            arg1 = str(self.end5File)
            arg2 = str(self.end3File)

            # upstream args
            arg3 = "{:.0f}".format(self.upstream[0]) if self.upstream else DUMMY_VALUE
            arg4 = "{:.0f}".format(self.upstream[1]) if self.upstream else DUMMY_VALUE
            sys.stdout.flush()
            retval = subprocess.call(["Rscript",
                                      ABSEQROOT + "/rscripts/masterScript.R",
                                      arg1,
                                      arg2,
                                      arg3,
                                      arg4],
                                     stdout=sys.stdout,
                                     stderr=sys.stdout
                                     )
            if retval != 0:
                print("-" * 30)
                print("Error detected in R plotting")
                print("-" * 30)
            # os.remove(PlotManager._tmpFile)        TODO: necessary?

    def processSingleInput(self, name, resultDirName=RESULT_FOLDER):
        """
        parses name into metadata as Rscript file for R plot to know where the same file is located.
        :param name: string
                    single sample name

        :param resultDirName: string
                    intermediate directory name (if any)
        :return: None
        """
        outDirName = os.path.join(self.args.outdir, resultDirName, name)
        self.metadata.append((outDirName, name))

    def requestedSamples(self, allFiles, resultDirName=RESULT_FOLDER, displayMatch=False):
        """
        return a filtered list of tuples according to allFiles and -rs argument; if -rs was specified, only samples
        requested will be returned with the name as provided within -rs, otherwise all samples in alLFiles are returned
        with AbSeq's inferred sample names

        :param allFiles: tuple
                        2-tuple of _R1 and _R2 files

        :param resultDirName: string
                        name of result folder (usually report/)

        :param displayMatch: bool
                        prints inferred pairing if set to true, to stdout

        :return: list of tuples
                        each tuple in the list has the format
                        [((sample1_r1, sample1_r2), sample_name1), (sample2, sample_name2), ... ]
                        where the first element has either a string (single ended) or a tuple (paired-end)
                        and the second element as the user provided sample name in -rs or AbSeq's inferred
                        name if -rs was not provided
        """
        from abseq.IgRepertoire.igRepUtils import inferSampleName
        # mapping = [ (sample1_r1, sample1_r2, rsName1), (sample2, rsName2), (...) ... ]
        mapping = []
        for sample in allFiles:
            if type(sample) == tuple:
                # paired end sample
                f1name, _ = sample
                merger = True
            else:
                f1name = sample
                merger = False

            canonicalName = inferSampleName(f1name, merger=merger, fastqc=(self.args.task.lower() == 'fastqc'))
            outDirName = os.path.join(self.args.outdir, resultDirName, canonicalName)
            self.metadata.append((outDirName, canonicalName))
            mapping.append((sample, canonicalName))

        if not self.rscriptArgs or self.rscriptArgs == 'off':
            return mapping

        # given directory names and canonical names, try to map files to names
        self._mapAllFiles()

        # metadata consists of old (abseq inferred) name, empty it and refill with user provided names
        self.metadata = []
        mapping = []
        #                  (to)
        # map inferred name => -rs name (user provided name)
        requestedSampleCanonicalNames = {}
        for rsName, val in self._nameToFileMap.items():
            requestedSampleCanonicalNames[val[-1]] = rsName

        if displayMatch:
            print("AbSeq inferred the following pairing from -rs argument:")
            for k, v in requestedSampleCanonicalNames.items():
                print("\t{:<20}: {:>20}".format(v, k))
            sys.stdout.flush()

        for sample in allFiles:
            if type(sample) == tuple:
                # paired end sample
                f1name, _ = sample
                merger = True
            else:
                f1name = sample
                merger = False
            canonicalName = inferSampleName(f1name, merger=merger, fastqc=(self.args.task.lower() == 'fastqc'))
            if canonicalName in requestedSampleCanonicalNames:
                userProvidedName = requestedSampleCanonicalNames[canonicalName]
                mapping.append((sample, userProvidedName))
                self.metadata.append((os.path.join(self.args.outdir, resultDirName, userProvidedName),
                                      userProvidedName))
        return mapping

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

    def _findBestMatch(self, sampleName, useProvidedName=False):
        if useProvidedName:
            tmp = list(self._nameToFileMap[sampleName])
            tmp[0] = os.path.join(os.path.split(tmp[0])[0], sampleName)
            tmp[1] = sampleName
            return tuple(tmp)
        return self._nameToFileMap[sampleName]

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

    def _mapAllFiles(self):
        if self.rscriptArgs and self.rscriptArgs != 'off':
            sampleNameMap = defaultdict(list)
            fileNameMap = defaultdict(list)
            namedSamples = list(set(itertools.chain(*self.rscriptArgs)))
            for rs in namedSamples:
                for dirName, sampleName in self.metadata:
                    sampleNameMap[rs].append(max(_nameMatch(dirName, rs), _nameMatch(dirName, rs)))
            for dirName, sampleName in self.metadata:
                for rs in namedSamples:
                    fileNameMap[sampleName].append(max(_nameMatch(dirName, rs), _nameMatch(dirName, rs)))
            taken = set()
            for rs in namedSamples:
                bestInd = int(np.argmax(sampleNameMap[rs]))
                # check scores of each (dirName, sampleName), if there's an obvious winner, we're done
                # also make sure that we've not chosen the winner before,
                # it doesn't make sense to have a one to many rel.s else, resolve ties
                if sampleNameMap[rs].count(sampleNameMap[rs][bestInd]) == 1 and self.metadata[bestInd][1] not in taken:
                    self._nameToFileMap[rs] = self.metadata[bestInd]
                    taken.add(self._nameToFileMap[rs][1])
                else:
                    # if there's a tie:
                    # 1. get sampleName of all tiedIndices
                    # 2. check that at least of of the (dir, sampleName) score has this sample as their max score
                    # 3. if it doesn't, raise ambiguous exception
                    # 4. if it does, make sure there's only one winner, then assign that as the match.

                    # 1.
                    tieInds = [i for i, x in enumerate(sampleNameMap[rs]) if x == sampleNameMap[rs][bestInd]]
                    assert len(tieInds) > 1 or self.metadata[bestInd][1] in taken
                    if self.metadata[bestInd][1] in taken:
                        raise Exception("Sample name {} conflicts with {}"
                                        .format(rs, self.metadata[bestInd][1]))
                    tiedSamples = [self.metadata[i][1] for i in tieInds if self.metadata[i][1] not in taken]
                    sampleIndex = namedSamples.index(rs)

                    # 2.
                    # check that either one of the files have MAX score = sample
                    for s in tiedSamples:
                        if sampleIndex == int(np.argmax(fileNameMap[s])):
                            break
                    else:  # 3
                        other = namedSamples[int(np.argmax(fileNameMap[s]))]
                        raise Exception(
                            "Your naming scheme for {} is too ambiguous with {} and {}, try a more specific one"
                            .format(s, other, rs))

                    # 4
                    # check that there's no more tie. You can't have a tie in a tiebreaker algorithm ...
                    if len(set([fileNameMap[s][sampleIndex] for s in tiedSamples])) > 1:
                        raise Exception("Multiple files match this sample name!")

                    bestSampleScoreInd = int(np.argmax([fileNameMap[s][sampleIndex] for s in tiedSamples]))
                    tmp = [i for i, p in enumerate(self.metadata) if p[1] == tiedSamples[bestSampleScoreInd]]
                    assert len(tmp) == 1
                    self._nameToFileMap[rs] = self.metadata[tmp[0]]
                    taken.add(self._nameToFileMap[rs][1])


def _nameMatch(string1, string2, deletionPenalty=-3, insertionPenalty=-3, matchScore=5, mismatchScore=-3):
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
    for i in range(1, len(string2) + 1):
        for j in range(1, len(string1) + 1):
            matrix[i][j] = max(matrix[i][j - 1] + deletionPenalty, matrix[i - 1][j] + insertionPenalty,
                               matrix[i - 1][j - 1] + (
                                   matchScore if string1[j - 1] == string2[i - 1] else mismatchScore), 0)
            maxval = max(maxval, matrix[i][j])
    return maxval
