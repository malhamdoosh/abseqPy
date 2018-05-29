'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
from __future__ import division
import gc
import os
import logging
import sys
import inspect

from collections import Counter, defaultdict
from os.path import exists
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas.io.parsers import read_csv
from pandas.io.pytables import read_hdf
from numpy import Inf, random, isnan, logical_not

from abseqPy.IgMultiRepertoire.AbSeqWorker import AbSeqWorker
from abseqPy.IgRepAuxiliary.upstreamAuxiliary import plotUpstreamLenDist, extractUpstreamSeqs, \
    writeCountsCategoriesToFile, findUpstreamMotifs
from abseqPy.IgRepAuxiliary.primerAuxiliary import addPrimerData, generatePrimerPlots
from abseqPy.config import FASTQC, AUX_FOLDER, HDF_FOLDER, DEFAULT_TASK, DEFAULT_MERGER, DEFAULT_TOP_CLONE_VALUE
from abseqPy.IgRepertoire.igRepUtils import compressCountsGeneLevel, gunzip, fastq2fasta, mergeReads, \
    writeListToFile, writeSummary, compressSeqGeneLevel, compressSeqFamilyLevel, \
    createIfNot, safeOpen, detectFileFormat, countSeqs
from abseqPy.versionManager import writeParams
from abseqPy.logger import printto, setupLogger, LEVEL
from abseqPy.IgRepAuxiliary.productivityAuxiliary import refineClonesAnnotation
from abseqPy.IgRepReporting.igRepPlots import plotSeqLenDist, plotSeqLenDistClasses, plotVenn, plotDist, eitherExists
from abseqPy.IgRepAuxiliary.annotateAuxiliary import annotateIGSeqRead
from abseqPy.IgRepReporting.abundanceReport import writeAbundanceToFiles
from abseqPy.IgRepReporting.productivityReport import generateProductivityReport
from abseqPy.IgRepReporting.diversityReport import generateDiversityReport, writeClonotypeDiversityRegionAnalysis
from abseqPy.IgRepAuxiliary.diversityAuxiliary import annotateSpectratypes, \
    annotateClonotypes
from abseqPy.IgRepAuxiliary.restrictionAuxiliary import findHits, \
    findHitsRegion, scanRestrictionSitesSimple, loadRestrictionSites
from abseqPy.IgRepReporting.restrictionReport import generateOverlapFigures
from abseqPy.utilities import hasLargeMem


# the following are conditionally imported in functions that require them to reduce abseq's dependency list
# It's here for a simple glance of required dependencies (generateMotifs uses TAMO)
# from IgRepAuxiliary.SeqUtils import generateMotifs
# from TAMO.Clustering.UPGMA import UPGMA
# from TAMO.Clustering.UPGMA import DFUNC
# from TAMO.Clustering.UPGMA import print_tree
# from TAMO.Clustering.UPGMA import create_tree_phylip
# from TAMO.Clustering.UPGMA import print_tree_id
# from TAMO import MotifTools


class IgRepertoire:
    """
    creates an AbSeq.IgRepertoire object with QC methods
    """
    def __init__(self, f1, f2=None, name=None, fmt=None, chain='hv', seqtype='dna', domainSystem='imgt',
                 merger=DEFAULT_MERGER, outdir='.', threads=1, bitscore=(0, Inf), alignlen=(0, Inf),
                 sstart=(1, Inf), qstart=(1, Inf), clonelimit=DEFAULT_TOP_CLONE_VALUE, detailedComposition=False,
                 actualqstart=-1, trim5=0, trim3=0,
                 fr4cut=True, primer=None, primer5endoffset=0, primer5end=None, primer3end=None,
                 upstream=None, sites=None, database="$IGBLASTDB", report_interim=False, task=DEFAULT_TASK, log=None,
                 yaml=None):
        """

        :param f1: string
                                path to read 1 file
        :param f2: string
                                path to read 2 file. This is optional
        :param name: string
                                name to refer this sample as
        :param fmt: string
                                accepted values are fasta, fa, fastq, fq. f1 and f2(if present) should have the
                                same format
        :param chain: string
                                accepted values are lv, hv, kv for lambda variable, heavy variable and kappa variable
                                respectively
        :param seqtype: string
                                accepted values are dna or protein
        :param domainSystem: string
                                imgt or kabat numbering
        :param merger: string
                                name of merger to use. This is ignored if f2 is not provided
        :param outdir: string

                                path to results directory. implicitly create if doesn't exist
        :param threads: string
                                number of threads to run this sample with
        :param bitscore: list or tuple
                                iterable and indexable of length 2 denoting the min and max value to use for
                                filtering sequences that do not fall within the provided range.
                                The bitscore filter applies to the V germline alignment only.
        :param alignlen: list or tuple
                                iterable and index-able of length 2 denoting the min and max value to use for
                                filtering sequences that do not fal within the provided range.
                                The alignlen filter applies to the V germline only
        :param sstart: list or tuple
                                iterable and index-able of length 2 denoting the min and max value to use for
                                filtering sequences that do not fall within the provided range.
                                The sstart filter applies to the V germline only. In this case, subject start
                                denotes the starting index of V germline when aligned to the query sequence
        :param qstart: list or tuple
                                iterable and indexable of length 2 denoting the min and max value to use for
                                filtering sequences that do not fall within the provided range.
                                The qstart filter applies to the V germline only. In this case, query start
                                denotes the starting index of the query sequence when aligned to the V germline gene
        :param clonelimit: int
                                number of CDR3 clones to output into
                                diversity/<sample_name>_clonotypes_<clonelimit>_[over|under].csv.gz
                                This csv file contains CDR3 AA sequences with their counts. Also accepts
                                np.Inf to retain all clones
        :param detailedComposition: bool
                                should composition logos be plotted with IGV gene segregation?
        :param actualqstart: int
                                number of nucleotides to ignore at the beginning of the sequence before
                                V germline starts aligning. Leave this as -1 to let AbSeq automatically infer
                                from IgBLAST's alignment. This argument has no effect when aligning 5'
                                primer during primer specificity analysis
        :param trim5: int
                                number of nucleotides to trim on the 5' end of V domain
                                This argument has no effect when aligning 5' primer during primer specificity analysis
        :param trim3: int or list of strings
                                number of nucleotides to trim on the 3' end of V domain
                                This argument has no effect when aligning 3' primer during primer specificity analysis.
                                If a list of strings was provided instead,
                                then the end of the sequences will be trimmed at the starting position (incl) of
                                one of the (best matched) sequence in trim3
        :param fr4cut: bool
                                fr4cut automatically cut sequence after end of J germline gene
                                (extend 3' end of J gene to get actual FR4 end position if mismatches occur). If this is
                                set to False, trimming of the 3' end will depend on trim3's option
        :param primer: int
                                (not implemented yet)
        :param primer5endoffset: int
                                number of nucleotides to offset before staring to align the 5' primer sequences. Only
                                used during primer specificity analysis
        :param primer5end: string
                                path to 5' primer FASTA file. Only required if task was primer
        :param primer3end: string
                                path to 3' primer FASTA file. Only required if task was primer
        :param upstream: list or tuple
                                iterable or index-able of length 2 that denotes the start and end position of upstream
                                sub-sequences.
        :param sites: string
                                path to restriction sites file. Required only if task was rsa or rsasimple
        :param database: string
                                path to IgBLAST database (directory should contain output of makeblastdb).
                                Environment variables are also accepted, for example, export IGBLASTDB=/path/to/db
                                will require db to be the string "$IGBLASTDB"
        :param report_interim: bool
                                create intermediate report (not implemented yet)
        :param task: string
                                all, annotate, abundance, diversity, fastqc, productivity,
                                primer, 5utr, rsasimple, rsa, seqlen, secretion, seqlenclass. This variable
                                is responsible for the "banner" printed in the log file.
        :param log: string
                                path to logger file
        :param yaml: string
                                dummy variable. Used in commandline mode
        """
        fargs, _, _, values = inspect.getargvalues(inspect.currentframe())
        self.args = dict([(arg, values[arg]) for arg in fargs if arg != 'self'])
        # todo
        # sanitizeArgs(self.args)

        self.task = task.lower().strip()
        self.chain = chain
        self.name = name
        self.fr4cut = fr4cut
        self.reportInterim = report_interim

        # directory creation
        outputDir = os.path.abspath(outdir)
        self.hdfDir = os.path.join(outputDir, HDF_FOLDER, self.name) + os.path.sep
        self.auxDir = os.path.join(outputDir, AUX_FOLDER, self.name) + os.path.sep

        if not os.path.exists(self.hdfDir):
            os.makedirs(self.hdfDir)
        if not os.path.exists(self.auxDir):
            os.makedirs(self.auxDir)

        self.threads = threads
        self.primer = primer
        self.db = os.path.abspath(os.path.expandvars(database))
        self.bitScore = bitscore
        self.clonelimit = clonelimit
        self.alignLen = alignlen
        self.sStart = sstart
        self.qStart = qstart
        self.seqType = seqtype
        self.domainSystem = domainSystem

        if task in ['secretion', '5utr']:
            self.upstream = upstream

        if task in ['rsa', 'rsasimple']:
            self.sitesFile = sites

        if task in ['diversity', 'all']:
            self.detailedComposition = detailedComposition

        self.actualQstart = actualqstart

        self.trim5End = trim5
        self.trim3End = trim3

        self.end5 = primer5end
        self.end3 = primer3end
        self.end5offset = primer5endoffset

        self.readFile1 = f1
        self.readFile2 = f2
        self.format = fmt if fmt is not None else detectFileFormat(self.readFile1)
        self.merger = merger
        self.merge = 'no' if self.merger is None else 'yes'

        self.seqsPerFile = int(10.0 ** 5 / 2)
        self.cloneAnnot = None
        self.cloneSeqs = None
        self.readFile = None
        # True of any of the following directories are already created. We need to distinguish this
        # from the beginning because AbSeq also re-reads HDF within the same analysis to prevent
        # pickling self.cloneAnnot, self.cloneSeqs into multiprocessing.Queue
        self.warnOldDir = any(map(lambda x: exists(os.path.join(self.hdfDir, x)),
                                  ["abundance", "productivity", "diversity", "restriction_sites",
                                   "primer_specificity", 'utr5', 'secretion']))

        setupLogger(self.name, self.task, log)
        writeParams(self.args, self.auxDir)
        self._tasks = []
        self._setupTasks()
        self._summaryFile = os.path.join(self.auxDir, "summary.txt")

    def runFastqc(self):
        logger = logging.getLogger(self.name)

        if self.format == 'fasta':
            printto(logger, "Fasta file extension detected, will not perform fastqc", LEVEL.WARN)
            return

        outDir = os.path.join(self.auxDir, "fastqc")

        if not os.path.isdir(outDir):
            os.makedirs(outDir)

        filename = os.path.join(outDir, self.readFile1.split(os.path.sep)[-1]
                                .replace(".fastq", "")
                                .replace(".gz", "") + "_fastqc.html")

        if os.path.exists(filename):
            printto(logger, "fastqc was already performed on this library.", LEVEL.WARN)
            return

        printto(logger, "Fastqc is running ... ")

        command = "%s -o %s -t %d %s"
        # check for presence of file2 before concatenating str and None(None when self.readFile2 is empty/not provided)
        os.system(command % (FASTQC, outDir, self.threads,
                             self.readFile1 + " " + (self.readFile2 if self.readFile2 is not None else "")))
        paramFile = writeParams(self.args, outDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)
        printto(logger, "Fastqc has completed.")

    def mergePairedReads(self):
        logger = logging.getLogger(self.name)
        if self.merge != 'yes':
            self.readFile = self.readFile1
        else:
            mergedFastq = mergeReads(self.readFile1, self.readFile2,
                                     self.threads, self.merger, self.hdfDir, stream=logger)
            self.readFile = mergedFastq

    def annotateClones(self, filterOutDir=None):
        """
        annotate clones from self.read using IgBLAST. self.cloneAnnot will be a dataframe
        with annotated clones
        :param filterOutDir: string.
                filtered clones will be placed in this directory under the
                name /filterOutDir/<self.name>_filtered_out_clones.txt
        :return: None
        """
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, "annot")
        outAuxDir = os.path.join(self.hdfDir, "annot")

        if not os.path.isdir(outResDir):
            os.makedirs(outResDir)
        if not os.path.isdir(outAuxDir):
            os.makedirs(outAuxDir)

        cloneAnnotFile = os.path.join(outAuxDir, self.name + "_clones_annot.h5")

        if self.readFile is None:
            self.mergePairedReads()

        writeSummary(self._summaryFile, "RawReads", countSeqs(self.readFile1))

        if exists(cloneAnnotFile):
            if self.task == "annotate":
                printto(logger, "\tClones annotation file found and no further work needed ... " +
                        os.path.basename(cloneAnnotFile))
            else:
                printto(logger, "\tClones annotation file found and being loaded ... " +
                        os.path.basename(cloneAnnotFile))
                self.cloneAnnot = read_hdf(cloneAnnotFile, "cloneAnnot")
            # write (update) number of annotated reads - this is possible because we always
            # save the unfiltered cloneannot dataframe, and re-filter after reloading
            writeSummary(self._summaryFile, "AnnotatedReads", self.cloneAnnot.shape[0])
        else:
            if not exists(self.readFile):
                raise Exception(self.readFile + " does not exist!")

            # Convert FASTQ file into FASTA format
            if self.format == 'fastq':
                readFasta = fastq2fasta(self.readFile, self.hdfDir, stream=logger)
            elif self.format == 'fasta':
                # unzip the fasta file if need be
                readFasta = gunzip(self.readFile)
            else:
                raise Exception('unknown file format! ' + self.format)
            #             if self.trim3End > 0 or self.trim5End > 0:
            #                 trimSequences(readFasta)
            #                 self.trimmed = True

            # Estimate the IGV family abundance for each library
            (self.cloneAnnot, filteredIDs) = annotateIGSeqRead(self, readFasta,
                                                               self.seqType, outdir=outAuxDir,
                                                               domainSystem=self.domainSystem,
                                                               stream=logger)
            sys.stdout.flush()
            gc.collect()

            if len(filteredIDs):
                writeListToFile(filteredIDs, os.path.join(outAuxDir, self.name + "_unmapped_clones.txt"))
            # export the CDR/FR annotation to a file
            printto(logger, "\tClones annotation file is being written to " +
                    os.path.basename(cloneAnnotFile))
            self.cloneAnnot.to_hdf(cloneAnnotFile, "cloneAnnot", mode='w')
            paramFile = writeParams(self.args, outResDir)
            printto(logger, "The analysis parameters have been written to " + paramFile)

            # write number of annotated reads
            writeSummary(self._summaryFile, "AnnotatedReads", self.cloneAnnot.shape[0])

        printto(logger, "Number of clones that are annotated is {0:,}".format(
                int(self.cloneAnnot.shape[0])), LEVEL.INFO)

        filterOutDir = outAuxDir if filterOutDir is None else filterOutDir
        # Filter clones based on bitscore, alignLen, qStart, and sStart
        selectedRows = self._cloneAnnotFilteredRows(logger)
        filteredIDs = self.cloneAnnot[logical_not(selectedRows)]

        if len(filteredIDs) > 0:
            filteredIDs = filteredIDs[['vgene', 'vstart', 'vqstart', 'bitscore', 'alignlen']]
            filteredIDs.to_csv(os.path.join(filterOutDir, self.name + "_filtered_out_clones.txt"),
                               sep="\t", header=True, index=True)

        retained = int(self.cloneAnnot.shape[0]) - len(filteredIDs)

        printto(logger, 'Percentage of retained clones is {:.2%} ({:,}/{:,})'.format(
                retained / self.cloneAnnot.shape[0],
                retained,
                int(self.cloneAnnot.shape[0])), LEVEL.INFO)

        self.cloneAnnot = self.cloneAnnot[selectedRows]

        # generate plot of clone sequence length distribution
        outputFile = os.path.join(outResDir, self.name + '_all_clones_len_dist.png')
        noOutlierOutputFile = os.path.join(outResDir, self.name + '_all_clones_len_dist_no_outliers.png')
        seqLengths = defaultdict(int)
        if not eitherExists(outputFile) or not eitherExists(noOutlierOutputFile):
            for record in SeqIO.parse(gunzip(self.readFile), self.format):
                if record.id in self.cloneAnnot.index:
                    seqLengths[len(record)] += 1

            count = Counter(seqLengths)

            plotSeqLenDist(count, self.name, outputFile, self.format,
                           maxbins=40, histtype='bar', removeOutliers=False,
                           normed=True, stream=logger)
            # generate plot of clone sequence length distribution with outliers removed
            plotSeqLenDist(count, self.name, noOutlierOutputFile, self.format,
                           maxbins=40, histtype='bar', removeOutliers=True,
                           normed=True, stream=logger)
        else:
            printto(logger, "File found ... {}".format(os.path.basename(outputFile)), LEVEL.WARN)
            printto(logger, "File found ... {}".format(os.path.basename(noOutlierOutputFile)), LEVEL.WARN)

        # finally, write number of filtered reads
        writeSummary(self._summaryFile, "FilteredReads", self.cloneAnnot.shape[0])

    def analyzeAbundance(self):
        # Estimate the IGV family abundance for each library
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, "abundance")
        outAuxDir = os.path.join(self.hdfDir, "abundance")

        createIfNot(outResDir)
        createIfNot(outAuxDir)

        if self.cloneAnnot is None:
            self.annotateClones(outAuxDir)

        writeAbundanceToFiles(self.cloneAnnot, self.name, outResDir, self.chain, stream=logger)
        gc.collect()
        paramFile = writeParams(self.args, outResDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzeProductivity(self, inplaceProductive=True, inplaceFiltered=True):
        """
        analyze sample productivity

        :param inplaceProductive:
                    if this is set to true, self.cloneAnnot and self.cloneSeqs will only contain
                    productive sequences after this method finishes.

        :param inplaceFiltered:
                    if this is set to true, self.cloneAnnot and self.cloneSeqs will only contain
                    unfiltered sequences after this method finishes.
        :return: None
        """
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, "productivity")
        outAuxDir = os.path.join(self.hdfDir, "productivity")

        createIfNot(outResDir)
        createIfNot(outAuxDir)

        refinedCloneAnnotFile = os.path.join(outAuxDir, self.name + "_refined_clones_annot.h5")
        cloneSeqFile = os.path.join(outAuxDir, self.name + "_clones_seq.h5")

        if not exists(refinedCloneAnnotFile) or not exists(cloneSeqFile):
            if self.cloneAnnot is None:
                self.annotateClones(outAuxDir)
            self._reloadAnnot()
            #             if self.trimmed:
            #                 self.trim3End = 0
            #                 self.trim5End = 0
            #             elif self.trim3End > 0 or self.trim5End > 0:
            #                 print("WARNING: if trimming was applied in the 'annotate' step"
            #                       ", you may not need trimming")
            # print(sys.getsizeof(self.cloneAnnot) / (1024.**3)) # in GB
            (self.cloneAnnot, self.cloneSeqs) = refineClonesAnnotation(outAuxDir, self.name,
                                                                       self.cloneAnnot, self.readFile,
                                                                       self.format, self.actualQstart,
                                                                       self.chain, self.fr4cut,
                                                                       self.trim5End, self.trim3End,
                                                                       self.seqsPerFile, self.threads, self.db,
                                                                       stream=logger)
            gc.collect()
            # if generateReport:
            # export the CDR/FR annotation to a file                
            printto(logger, "The refined clone annotation file is being written to "
                    + os.path.basename(refinedCloneAnnotFile))
            self.cloneAnnot.to_hdf(refinedCloneAnnotFile, "refinedCloneAnnot", mode='w', complib='blosc')

            printto(logger, "The clone protein sequences are being written to " + os.path.basename(cloneSeqFile))
            self.cloneSeqs.to_hdf(cloneSeqFile, "cloneSequences", mode='w', complib='blosc')

            paramFile = writeParams(self.args, outResDir)
            printto(logger, "The analysis parameters have been written to " + paramFile)
            # although self.cloneAnnot is already filtered,
            # reapply filtering because vqstart might've changed post refinement
            printto(logger, "Applying filtering criteria to refined datafames")
        else:
            printto(logger, "The refined clone annotation files were found and being loaded ... " +
                    os.path.basename(refinedCloneAnnotFile))

            self.cloneAnnot = read_hdf(refinedCloneAnnotFile, "refinedCloneAnnot")
            printto(logger, "\tClone annotation was loaded successfully")

            self.cloneSeqs = read_hdf(cloneSeqFile, "cloneSequences")
            printto(logger, "\tClone sequences were loaded successfully")

            # since we loaded it from the saved (old) HDF5 dataframes, we need to re-apply all filtering criteria
            printto(logger, "\tApplying filtering criteria to loaded HDF5 dataframes")

        before = self.cloneAnnot.shape[0]
        selectedRows = self._cloneAnnotFilteredRows(logger)
        self.cloneAnnot = self.cloneAnnot[selectedRows]
        self.cloneSeqs = self.cloneSeqs.loc[self.cloneAnnot.index]
        printto(logger, "\tPercentage of retained clones is {:.2%} ({:,}/{:,})"
                .format(self.cloneAnnot.shape[0] / before, self.cloneAnnot.shape[0], before))

        # update number of filtered reads "POST REFINEMENT"
        writeSummary(self._summaryFile, "FilteredReads", self.cloneAnnot.shape[0])

        # display statistics
        printto(logger, "Productivity report is being generated ... ")
        generateProductivityReport(self.cloneAnnot, self.cloneSeqs, self.name, self.chain, outResDir, stream=logger)

        before = int(self.cloneAnnot.shape[0])
        inFrame = self.cloneAnnot[self.cloneAnnot['v-jframe'] == 'In-frame']

        # do not filter out "filtered" yet! - that has nothing to do with productivity
        if inplaceProductive:
            cloneAnnot = self.cloneAnnot = inFrame[inFrame['stopcodon'] == 'No']
            self.cloneSeqs = self.cloneSeqs.loc[self.cloneAnnot.index]
        else:
            cloneAnnot = inFrame[inFrame['stopcodon'] == 'No']
        printto(logger, "Percentage of productive clones {:.2%} ({:,}/{:,})".format(
            0 if before == 0 else cloneAnnot.shape[0] / before,
            int(cloneAnnot.shape[0]),
            int(before)
            ), LEVEL.INFO)

        # write number of productive reads
        writeSummary(self._summaryFile, "ProductiveReads", cloneAnnot.shape[0])

        # filter out "filtered" now
        if inplaceFiltered:
            self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['filtered'] == 'No']
            self.cloneSeqs = self.cloneSeqs.loc[self.cloneAnnot.index]

    def analyzeDiversity(self):
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, "diversity")
        outAuxDir = os.path.join(self.hdfDir, "diversity")

        if self.cloneAnnot is None or self.cloneSeqs is None:
            # we analyze productive clones ONLY
            self.analyzeProductivity(inplaceProductive=True, inplaceFiltered=True)

        if len(self.cloneAnnot) == 0:
            printto(logger, "WARNING: There are no productive sequences found (post-refinement) in {},"
                            " skipping diversity analysis.".format(self.name), LEVEL.WARN)
            return

        createIfNot(outResDir)
        createIfNot(outAuxDir)

        gc.collect()

        # Identify spectratypes 
        printto(logger, "Spectratypes are being calculated ... ")
        spectraTypes = annotateSpectratypes(self.cloneAnnot, amino=True)

        # Identify clonotypes 
        printto(logger, "Clonotypes are being generated ... ")
        clonoTypes = annotateClonotypes(self.cloneSeqs, segregate=self.detailedComposition, removeNone=True)

        generateDiversityReport(spectraTypes, clonoTypes, self.name, outResDir, self.clonelimit,
                                threads=self.threads, segregate=self.detailedComposition, stream=logger)

        # todo: remove this for now - it's unoptimized and extremely slow
        # writeClonotypeDiversityRegionAnalysis(self.cloneSeqs, self.name, outResDir, stream=logger)

        paramFile = writeParams(self.args, outResDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzeRestrictionSitesSimple(self):
        # TODO: parallelize this function to run faster
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, "restriction_sites")
        outAuxDir = os.path.join(self.hdfDir, "restriction_sites")

        if not os.path.isdir(outResDir):
            os.makedirs(outResDir)

        if not os.path.isdir(outAuxDir):
            os.makedirs(outAuxDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: remove the 'restriction_sites' directory if you changed the filtering criteria.",
                    LEVEL.WARN)

        siteHitsFile = os.path.join(outResDir, self.name + "_{}_rsasimple.csv"
                                    .format(os.path.splitext(os.path.basename(self.sitesFile))[0]))
        overlap2File = siteHitsFile.replace('.csv', '_overlap_order2.csv')

        if exists(siteHitsFile):
            printto(logger, "Restriction sites were already scanned at ... " +
                    os.path.basename(siteHitsFile), LEVEL.WARN)
            rsaResults = read_csv(siteHitsFile, header=0)
            if exists(overlap2File):
                overlapResults = {}
                overlapResults['order2'] = read_csv(overlap2File, header=0, index_col=0)
            else:
                overlapResults = None
        else:
            self.annotateClones(outAuxDir)
            (rsaResults, overlapResults) = scanRestrictionSitesSimple(self.name,
                                                                      self.readFile, self.format,
                                                                      self.cloneAnnot, self.sitesFile,
                                                                      self.threads)
            rsaResults.to_csv(siteHitsFile,
                              header=True,
                              index=False)
            printto(logger, "RSA results were written to " + os.path.basename(siteHitsFile))
            if overlapResults.get("order2", None) is not None:
                overlapResults["order2"].to_csv(overlap2File,
                                                header=True, index=True)
        # # print out the results        
        generateOverlapFigures(overlapResults,
                               rsaResults.loc[rsaResults.shape[0] - 1, "No.Molecules"],
                               self.name, siteHitsFile, stream=logger)
        paramFile = writeParams(self.args, outResDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzeRestrictionSites(self):
        # todo
        raise NotImplementedError
        # logger = logging.getLogger(self.name)
        #
        # outResDir = os.path.join(self.auxDir, "restriction_sites")
        # outAuxDir = os.path.join(self.hdfDir, "restriction_sites")
        #
        # if not os.path.isdir(outResDir):
        #     os.makedirs(outResDir)
        #
        # if not os.path.isdir(outAuxDir):
        #     os.makedirs(outAuxDir)
        # elif self.warnOldDir:
        #     printto(logger, "WARNING: remove the 'restriction_sites' directory if you changed the filtering criteria.",
        #             LEVEL.WARN)
        #
        # siteHitsFile = os.path.join(outResDir, self.name + "_{}.csv"
        #                             .format(os.path.splitext(os.path.basename(self.sitesFile))[0]))
        #
        # if exists(siteHitsFile):
        #     print("Restriction sites were already searched at ... " + os.path.basename(siteHitsFile))
        #     return
        #
        # if self.cloneAnnot is None or self.cloneSeqs is None:
        #     self.analyzeProductivity(inplaceProductive=True, inplaceFiltered=False)
        #
        # rsites = loadRestrictionSites(self.sitesFile, stream=logger)
        # printto(logger, "Restriction sites are being searched ... ")
        # gc.collect()
        #
        # queryIds = self.cloneAnnot.index
        # siteHitsCount = {}
        # siteHitSeqsCount = {}
        # hitRegion = {}
        # siteHitSeqsGermline = {}
        # seqsCutByAny = 0
        # siteHitsSeqsIDs = {}
        # siteHitsSeqsIGV = {}
        # for site in rsites.keys():
        #     siteHitsCount[site] = 0
        #     siteHitSeqsCount[site] = 0
        #     hitRegion[site] = Counter({'fr1': 0, 'cdr1': 0,
        #                                'fr2': 0, 'cdr2': 0,
        #                                'fr3': 0, 'cdr3': 0,
        #                                'fr4': 0})
        #     siteHitSeqsGermline[site] = []
        #     siteHitsSeqsIDs[site] = set()
        #     siteHitsSeqsIGV[site] = set()
        # germline = {'fr1', 'fr2', 'fr3', 'cdr1', 'cdr2'}
        # procSeqs = 0
        # #         if hasLargeMem():
        # #             TODO: remember to make sure SeqIO.parse is parsing a unzipped self.readFile1
        # #                   (use safeOpen from IgRepertoire.utils) if not sure
        # #             records = SeqIO.to_dict(SeqIO.parse(self.readFile1, self.format))
        # #         else:
        # # SeqIO.index can only open string file names and they must be uncompressed
        # records = SeqIO.index(gunzip(self.readFile1), self.format)
        # for id_ in queryIds:
        #     record = records[id_]
        #     try:
        #         qsRec = self.cloneAnnot.loc[record.id].to_dict()
        #         qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
        #         if isnan(qsRec['fr4.end']):
        #             end = len(record.seq)
        #         else:
        #             end = int(qsRec['fr4.end'])
        #         seq = str(record.seq[qstart:end])
        #         seqRC = str(Seq(seq).reverse_complement())
        #         cut = False
        #         for site in siteHitsCount.keys():
        #             hits = findHits(seq, rsites[site])
        #             strand = "forward"
        #             if len(hits) == 0:
        #                 hits = findHits(seqRC, rsites[site])
        #                 strand = "reversed"
        #             if len(hits) > 0:
        #                 siteHitsCount[site] += len(hits)
        #                 siteHitSeqsCount[site] += 1
        #                 hitsRegion = findHitsRegion(qsRec, hits)
        #                 if (len(set(hitsRegion).intersection(germline)) > 0
        #                         and len(siteHitSeqsGermline[site]) < 10000):
        #                     siteHitSeqsGermline[site].append((strand, str(record.seq)))
        #                     siteHitsSeqsIGV[site].add(qsRec['vgene'].split('*')[0])
        #                 hitRegion[site] += Counter(hitsRegion)
        #                 siteHitsSeqsIDs[site].add(record.id)
        #                 cut = True
        #         if cut:
        #             seqsCutByAny += 1
        #         procSeqs += 1
        #         if procSeqs % self.seqsPerFile == 0:
        #             print('{}/{} sequences have been searched ... '.format(procSeqs, len(queryIds)))
        #     #                 break
        #     except BaseException as e:
        #         print(qstart, end, len(record.seq), str(record.seq))
        #         print(e)
        #         raise
        # print('{}/{} sequences have been searched ... '.format(procSeqs, len(queryIds)))
        # # # print out the results
        # f = open(siteHitsFile, 'w')
        # f.write("Enzyme,Restriction Site,No.Hits,Percentage of Hits (%),"
        #         "No.Molecules,Percentage of Molecules (%),FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4, V Germlines \n")
        # sites = sorted(siteHitSeqsCount, key=siteHitSeqsCount.get)
        # for site in sites:
        #     f.write("{},{},{},{:.3%},{},{:.3%},{},{},{},{},{},{},{},{}\n"
        #             .format(site, rsites[site], siteHitsCount[site],
        #                     siteHitsCount[site] / sum(siteHitsCount.values()),
        #                     siteHitSeqsCount[site], siteHitSeqsCount[site] / len(queryIds),
        #                     hitRegion[site]['fr1'], hitRegion[site]['cdr1'], hitRegion[site]['fr2'],
        #                     hitRegion[site]['cdr2'], hitRegion[site]['fr3'], hitRegion[site]['cdr3'],
        #                     hitRegion[site]['fr4'], '|'.join(siteHitsSeqsIGV[site])))
        #     # write the first 100 sequences cut in the germline of each restriction enzyme
        #     seqs = []
        #     for (strand, seq) in siteHitSeqsGermline[site]:
        #         seqs.append(SeqRecord(Seq(seq), id='seq' + str(len(seqs)) + strand))
        #     SeqIO.write(seqs, siteHitsFile.replace('.csv', '_germline' + site + '.fasta'), 'fasta')
        # f.write("Sequences cut by any of the above enzymes, {}, {:.3%}\n"
        #         .format(seqsCutByAny, seqsCutByAny / len(queryIds)))
        # f.close()
        # # Ven Diagram of overlapping sequences
        # plotVenn(siteHitsSeqsIDs, siteHitsFile.replace('.csv', '_venn.png'), stream=logger)
        # print("Restriction enzyme results were written to " + siteHitsFile)

    def analyzeSecretionSignal(self):
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, 'secretion')
        outAuxDir = os.path.join(self.hdfDir, 'secretion')

        if not os.path.exists(outResDir):
            os.makedirs(outResDir)

        if not os.path.exists(outAuxDir):
            os.makedirs(outAuxDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: Remove 'secretion' directory if you've changed the filtering criteria.",
                    LEVEL.WARN)

        # need self.cloneAnnot dataframe for further analysis
        if self.cloneAnnot is None:
            self.annotateClones(outAuxDir)

        printto(logger, "The diversity of the upstream of IGV genes is being analyzed ... ")

        upstreamFile = os.path.join(outAuxDir, self.name + "_secsig_{:.0f}_{:.0f}.fasta"\
                                    .format(self.upstream[0], self.upstream[1]))

        if not exists(upstreamFile):
            extractUpstreamSeqs(self.cloneAnnot, self.readFile, self.upstream, upstreamFile, stream=logger)
        else:
            printto(logger, "\tUpstream sequences file was found! ... " + os.path.basename(upstreamFile), LEVEL.WARN)

        upstreamFile = os.path.abspath(upstreamFile)

        expectLength = self.upstream[1] - self.upstream[0] + 1

        # plot the distribution of sequence length
        plotUpstreamLenDist(upstreamFile, expectLength, self.name, outResDir, stream=logger)

        if expectLength != Inf:
            # classify secretion signals based on length, ATG location, gene and gene family

            # ----------------------------------------------------------------
            #                   analyze intact secretion signals
            # ----------------------------------------------------------------
            #  this means expectLength[0] == expectLength[1] (sequences with exactly expectLength in length only)
            printto(logger, "\tAnalyzing intact secretion signals", LEVEL.DEBUG)
            for level in ['variant', 'gene', 'family']:
                findUpstreamMotifs(upstreamFile, self.name, outAuxDir, outResDir,
                                   [expectLength, expectLength], level=level, startCodon=True,
                                   threads=self.threads, stream=logger)

            # ----------------------------------------------------------------
            #                 analyze trimmed secretion signals
            # ----------------------------------------------------------------
            printto(logger, "\tAnalyzing trimmed secretion signals", LEVEL.DEBUG)
            for level in ['variant', 'gene', 'family']:
                findUpstreamMotifs(upstreamFile, self.name, outAuxDir, outResDir, [1, expectLength - 1], level=level,
                                   startCodon=True, threads=self.threads, stream=logger)

        paramFile = writeParams(self.args, outResDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyze5UTR(self):
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, 'utr5')
        outAuxDir = os.path.join(self.hdfDir, 'utr5')

        if not os.path.exists(outResDir):
            os.makedirs(outResDir)

        if not os.path.exists(outAuxDir):
            os.makedirs(outAuxDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: Remove 'utr5' directory if you've changed the filtering criteria",
                    LEVEL.WARN)

        # requires self.cloneAnnot dataframe for further analysis
        if self.cloneAnnot is None:
            self.annotateClones(outAuxDir)

        printto(logger, "The diversity of the upstream of IGV genes is being analyzed ... ")

        upstreamFile = os.path.join(outAuxDir, self.name + "_5utr_{:.0f}_{:.0f}.fasta"\
                                    .format(self.upstream[0], self.upstream[1]))

        if not exists(upstreamFile):
            extractUpstreamSeqs(self.cloneAnnot, self.readFile, self.upstream, upstreamFile, stream=logger)
        else:
            printto(logger, "\tUpstream sequences file was found! ... " + os.path.basename(upstreamFile), LEVEL.WARN)

        upstreamFile = os.path.abspath(upstreamFile)
        expectLength = self.upstream[1] - self.upstream[0] + 1

        plotUpstreamLenDist(upstreamFile, expectLength, self.name, outResDir, stream=logger)

        # if user provided values to upstream (and it's not Inf)
        if expectLength != Inf:
            # ----------------------------------------------------------------
            #                 analyze intact secretion signals
            # ----------------------------------------------------------------
            #  this means expectLength[0] == expectLength[1] (sequences with exactly expectLength in length only)
            for level in ['variant', 'gene', 'family']:
                findUpstreamMotifs(upstreamFile, self.name, outAuxDir, outResDir, [expectLength, expectLength],
                                   level=level, startCodon=True, type='5utr', clusterMotifs=True,
                                   threads=self.threads, stream=logger)

        paramFile = writeParams(self.args, outResDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzePrimerSpecificity(self):
        logger = logging.getLogger(self.name)

        outResDir = os.path.join(self.auxDir, 'primer_specificity')
        outAuxDir = os.path.join(self.hdfDir, 'primer_specificity')

        createIfNot(outResDir)
        createIfNot(outAuxDir)

        primerAnnotFile = os.path.join(outAuxDir, self.name + "_primer_annot.h5")

        # if we can't find hdf file, create it, else read it
        if not exists(primerAnnotFile):
            # Load self.cloneAnnot for further analysis.
            # skip checking for existence of dataframes, analyzeProd/Abun will do it for us
            if self.cloneAnnot is None:
                self.annotateClones(outAuxDir)

            # addPrimerData on an unfiltered self.cloneAnnot
            self._reloadAnnot()

            # add additional primer related data to the dataframe generated by either abundance/productivity analysis
            # before we begin primer analysis
            self.cloneAnnot = addPrimerData(self.cloneAnnot, self.readFile, self.format, self.fr4cut,
                                            self.trim5End, self.trim3End, self.actualQstart,
                                            self.end5, self.end3, self.end5offset, self.threads, stream=logger)
            # save new "primer column-ed dataframe" into primer_specificity directory
            self.cloneAnnot.to_hdf(primerAnnotFile, "primerCloneAnnot", mode='w', complib='blosc')

            # now we can safely apply the filter on self.cloneAnnot
            printto(logger, "\tApplying filtering criteria to primer specificity analysis dataframes")
        else:
            printto(logger, "The primer clone annotation files were found and being loaded ... ", LEVEL.WARN)
            self.cloneAnnot = read_hdf(primerAnnotFile, "primerCloneAnnot")
            printto(logger, "\tPrimer clone annotation loaded successfully")

            # since we loaded it from the saved (old) HDF5 dataframes, we need to re-apply all filtering criteria
            printto(logger, "\tApplying filtering criteria to loaded HDF5 dataframes")

        before = self.cloneAnnot.shape[0]
        selectedRows = self._cloneAnnotFilteredRows(logger)
        self.cloneAnnot = self.cloneAnnot[selectedRows]
        self.cloneSeqs = self.cloneSeqs.loc[self.cloneAnnot.index]
        printto(logger, "\tPercentage of retained clones is {:.2%} ({:,}/{:,})"
                .format(self.cloneAnnot.shape[0] / before, self.cloneAnnot.shape[0], before))

        # TODO: Fri Feb 23 17:13:09 AEDT 2018
        # TODO: check findBestMatchAlignment of primer specificity best match, see if align.localxx is used correctly!
        generatePrimerPlots(self.cloneAnnot, outResDir, self.name, self.end5, self.end3, stream=logger)

        paramFile = writeParams(self.args, outResDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzeSeqLen(self, klass=False):
        logger = logging.getLogger(self.name)

        outResdir = os.path.join(self.auxDir, 'annot')

        printto(logger, "Sequence {}length distribution is being calculated ... ".format('class ' if klass else ''))

        if not os.path.exists(outResdir):
            os.makedirs(outResdir)

        if klass:
            outputFile = os.path.join(outResdir, self.name + '_length_dist_classes.png')
            plotSeqLenDistClasses(self.readFile, self.name, outputFile, self.format, stream=logger)
        else:
            outputFile = os.path.join(outResdir, self.name + '_seq_length_dist.png')
            plotSeqLenDist(self.readFile, self.name, outputFile, self.format, maxbins=-1, stream=logger)

        paramFile = writeParams(self.args, outResdir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def _cloneAnnotFilteredRows(self, logger):
        printto(logger, "Clones are being filtered based on the following criteria: ", LEVEL.INFO)
        printto(logger, "\tBit score: " + repr(self.bitScore), LEVEL.INFO)
        printto(logger, "\tAlignment length: " + repr(self.alignLen), LEVEL.INFO)
        printto(logger, "\tSubject V gene start: " + repr(self.sStart), LEVEL.INFO)
        printto(logger, "\tQuery V gene start: " + repr(self.qStart), LEVEL.INFO)
        selectedRows = (
                (self.cloneAnnot['bitscore'] >= self.bitScore[0]) &     # check bit-Score
                (self.cloneAnnot['bitscore'] <= self.bitScore[1]) &
                (self.cloneAnnot['alignlen'] >= self.alignLen[0]) &     # check alignment length
                (self.cloneAnnot['alignlen'] <= self.alignLen[1]) &
                (self.cloneAnnot['vstart'] >= self.sStart[0]) &         # check subject (V gene) start position
                (self.cloneAnnot['vstart'] <= self.sStart[1]) &
                (self.cloneAnnot['vqstart'] >= self.qStart[0]) &        # check query (V gene) start position
                (self.cloneAnnot['vqstart'] <= self.qStart[1])
        )
        return selectedRows

    def _reloadAnnot(self):
        """
        on analysis that takes self.cloneAnnot as an argument and produces results based on it might require that
        the data frame be unfiltered (otherwise analysis conducted after it will all be on the filtered dataframe).
        For example, analyzeProductivity and analyzePrimerSpecificity uses self.cloneAnnot as an input and produces
        a refined dataframe from it. That refined dataframe should be unfiltered so that further analysis can
        customize filtering later on (and users can use the same HDF5 file for different filtering criteria)
        :return: None
        """
        aux = os.path.join(self.hdfDir, "annot")
        cloneAnnotFile = os.path.join(aux, self.name + "_clones_annot.h5")
        if os.path.exists(cloneAnnotFile):
            self.cloneAnnot = read_hdf(cloneAnnotFile, "cloneAnnot")
        else:
            raise Exception("Cannot reload self.cloneAnnot, file {} not found".format(cloneAnnotFile))

    def _minimize(self):
        # XXX: cloneAnnot and cloneSeqs are the largest objects in a IgReportoire object,
        # we remove them so that we can pickle them into the queue again.
        # When needed, these files will be loaded automatically later on anyway
        self.cloneAnnot = None
        self.cloneSeqs = None

    def _nextTask(self):
        if len(self._tasks) > 0:
            pack = self._tasks.pop()
            if type(pack) == str:
                return pack, [], {}
            else:
                # type(pack) = tuple: (str, dict) - see SeqLenClass
                return pack[0], [], pack[1]
        else:
            return None, [], {}

    def _setupTasks(self):
        logger = logging.getLogger(self.name)
        if self.task == 'all':
            self._tasks = [AbSeqWorker.FASTQC, AbSeqWorker.ANNOT, AbSeqWorker.ABUN, AbSeqWorker.PROD, AbSeqWorker.DIVER]
        elif self.task == 'fastqc':
            self._tasks = [AbSeqWorker.FASTQC]
        elif self.task == 'annotate':
            self._tasks = [AbSeqWorker.ANNOT]
        elif self.task == 'abundance':
            self._tasks = [AbSeqWorker.ABUN]
        elif self.task == 'productivity':
            self._tasks = [AbSeqWorker.PROD]
        elif self.task == 'diversity':
            self._tasks = [AbSeqWorker.DIVER]
        elif self.task == 'secretion':
            self._tasks = [AbSeqWorker.SECR]
        elif self.task == '5utr':
            self._tasks = [AbSeqWorker.UTR5]
        elif self.task == 'rsasimple':
            self._tasks = [AbSeqWorker.RSAS]
        elif self.task == 'rsa':
            self._tasks = [AbSeqWorker.RSA]
        elif self.task == 'primer':
            self._tasks = [AbSeqWorker.PRIM]
        elif self.task == 'seqlen':
            self._tasks = [AbSeqWorker.SEQLEN]
        elif self.task == 'seqlenclass':
            self._tasks = [(AbSeqWorker.SEQLEN, {'klass': True})]
        else:
            raise ValueError("Unknown task requested: {}".format(self.task))

        # make sure that if user specified either one of primer end file, we unconditionally run primer analysis
        # (duh)
        if self.task != 'primer' and (self.end3 or self.end5):
            printto(logger, "Primer file detected, conducting primer specificity analysis ... ", LEVEL.INFO)
            self._tasks.append(AbSeqWorker.PRIM)

        # easier to pop
        self._tasks = self._tasks[::-1]

    def analyzeIgProtein(self):
        raise NotImplementedError
        # sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'
        # sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        # self.cloneSeqs = read_csv(cloneSeqFile, sep='\t',
        #                                header=0, index_col=0)
        # self.readFile1 = self.outputDir + self.name
        # self.readFile1 += '_productive_prot.fasta'
        # if (not exists(self.readFile1)):
        #     print("Protein sequences are being prepared ...")
        #     records = []
        #     procSeqs = 0
        #     open(self.readFile1, 'w').close()
        #     for id in self.cloneSeqs.index:
        #         seq = ''.join(self.cloneSeqs.loc[id,].tolist()[1:])
        #         if '*' in seq:
        #             seq = seq.replace('*', 'X')
        #         rec = SeqRecord(Seq(seq), id=id, name="", description="")
        #         records.append(rec)
        #         procSeqs += 1
        #         if procSeqs % self.seqsPerFile == 0:
        #             print('\t{}/{} sequences have been processed ...  '.format(procSeqs, len(self.cloneSeqs.index)))
        #             sys.stdout.flush()
        #             SeqIO.write(records, open(self.readFile1, 'a'), 'fasta')
        #             records = []
        #     SeqIO.write(records, open(self.readFile1, 'a'), 'fasta')
        #     del records
        # else:
        #     print("File found ... " + os.path.basename(self.readFile1))
        # self.format = 'fasta'
        # self.readFile2 = None
        # self.seqType = 'protein'
        # self.bitScore = [0, Inf]
        # self.alignLen = [0, Inf]
        # self.sStart = [1, Inf]
        # if exists(self.outputDir + "/abundance/"):
        #     print("Protein sequences have been already analyzed ... ")
        # else:
        #     self.analyzeAbundance()


#     def extractProductiveRNAs(self):
# #         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'
# #         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
#         # v-j rearrangement frame distribution
#         vjframeDist = Counter(self.cloneAnnot['v-jframe'].tolist())
#         if NaN in vjframeDist.keys():
#             nanCounts = vjframeDist[NaN]
#             vjframeDist = Counter({'In-frame': vjframeDist['In-frame'],
#                                    'Out-of-frame': vjframeDist['Out-of-frame'] + nanCounts})
#         plotDist(vjframeDist, self.name, self.outputDir + self.name +
#                  '_vjframe_dist.png', title='V-D-J Rearrangement',
#                  proportion=False, rotateLabels=False)
#         print(vjframeDist)
#         del vjframeDist
#         if self.end5:
#             print("5-end analysis of all clones ... ")
#             self.write5EndPrimerStats(self.cloneAnnot, self.outputDir+self.name+
#                                       '_all_5end_')
#             invalid5Clones = self.cloneAnnot.index[self.cloneAnnot['5end'] == 'Indelled'].tolist()
#         if self.end3:
#             valid3End = Counter(self.cloneAnnot['3end'].tolist())
#             plotDist(valid3End, self.name, self.outputDir + self.name +
#                  '_all_3end_integrity_dist.png', title='Integrity of 3`-end Primer Sequence',
#                  proportion=True, rotateLabels=False)
#             invalid3Clones = self.cloneAnnot.index[self.cloneAnnot['3end'] == 'Indelled'].tolist()
#             print("Example of Indelled 3`-end:", invalid3Clones[1:10])
#             try:
#                 plotVenn({'5`-end':set(invalid5Clones), '3`-end':set(invalid3Clones)},
#                           self.outputDir + self.name +
#                      '_all_invalid_primers.png')
#                 del invalid5Clones, invalid3Clones
#             except:
#                 pass
#             del valid3End
#
#         OutOfFrame = self.cloneAnnot[self.cloneAnnot['v-jframe'] != 'In-frame']
#         OutOfFrameFamilyDist = compressCountsFamilyLevel(Counter(OutOfFrame['vgene'].tolist()))
#         plotDist(OutOfFrameFamilyDist, self.name, self.outputDir + self.name +
#                  '_notinframe_igv_dist.png',
#                   title='IGV Abundance of Not In-frame Sequences',
#                  proportion=True)
#         del OutOfFrameFamilyDist
#         OutOfFrame = OutOfFrame[OutOfFrame['v-jframe'] == 'Out-of-frame']
#         cdrLength = (OutOfFrame['cdr1.end'] - OutOfFrame['cdr1.start'] + 1) / 3
#         cdrLength = cdrLength.tolist()
#         histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name +
#                  '_outframe_cdr1_len_dist.png', dna=False,
#                   seqName='CDR1', normed=True, maxbins=10)
#         cdrGaps = Counter(OutOfFrame['cdr1.gaps'].tolist())
#         plotDist(cdrGaps, self.name, self.outputDir + self.name +
#                  '_outframe_cdr1_gaps_dist.png', title='Gaps in CDR1',
#                  proportion=False, rotateLabels=False)
#         frGaps = Counter(OutOfFrame['fr1.gaps'].tolist())
#         plotDist(frGaps, self.name, self.outputDir + self.name +
#                  '_outframe_fr1_gaps_dist.png', title='Gaps in FR1',
#                  proportion=False, rotateLabels=False)
#         del cdrLength, cdrGaps, frGaps
#         cdrLength = (OutOfFrame['cdr2.end'] - OutOfFrame['cdr2.start'] + 1) / 3
#         cdrLength = cdrLength.tolist()
#         histcals = plotSeqLenDist(cdrLength, self.name, self.outputDir + self.name +
#                  '_outframe_cdr2_len_dist.png', dna=False,
#                   seqName='CDR2', normed=True, maxbins=10)
#         cdrGaps = Counter(OutOfFrame['cdr2.gaps'].tolist())
#         plotDist(cdrGaps, self.name, self.outputDir + self.name +
#                  '_outframe_cdr2_gaps_dist.png', title='Gaps in CDR2',
#                  proportion=False, rotateLabels=False)
#         frGaps = Counter(OutOfFrame['fr2.gaps'].tolist())
#         plotDist(frGaps, self.name, self.outputDir + self.name +
#                  '_outframe_fr2_gaps_dist.png', title='Gaps in FR2',
#                  proportion=False, rotateLabels=False)
#         del cdrLength, cdrGaps, frGaps
#         cdrGaps = Counter([x if not isnan(x) else 'NA' for x in OutOfFrame['cdr3g.gaps'] ])
# #         print(len(cdrGaps))
#         plotDist(cdrGaps, self.name, self.outputDir + self.name +
#                  '_outframe_cdr3_gaps_dist.png', title='Gaps in Germline CDR3',
#                  proportion=False, rotateLabels=False)
#         frGaps = Counter(OutOfFrame['fr3g.gaps'].tolist())
#         plotDist(frGaps, self.name, self.outputDir + self.name +
#                  '_outframe_fr3_gaps_dist.png', title='Gaps in FR3 (Germline)',
#                  proportion=False, rotateLabels=False)
#         del cdrGaps, frGaps
#         if self.end5:
#             print("5-end analysis of out-of-frame clones ... ")
#             self.write5EndPrimerStats(OutOfFrame, self.outputDir+self.name+
#                                       '_outframe_5end_', 'Out-of-frame')
#             invalid5Clones = OutOfFrame.index[OutOfFrame['5end'] == 'Indelled'].tolist()
#         if self.end3:
#             valid3End = Counter(OutOfFrame['3end'].tolist())
#             plotDist(valid3End, self.name, self.outputDir + self.name +
#                  '_outframe_3end_integrity_dist.png', title='Integrity of 3`-end Primer Sequence(Out-of-frame)',
#                  proportion=True, rotateLabels=False)
#             invalid3Clones = OutOfFrame.index[OutOfFrame['3end'] == 'Indelled'].tolist()
#             print("Example of out-of-frame Indelled 3`-end:", invalid3Clones[1:10])
#             print("Example of out-of-frame valid 3`-end:", OutOfFrame.index[OutOfFrame['3end'] != 'Indelled'].tolist()[1:10])
#             try:
#                 plotVenn({'5`-end':set(invalid5Clones), '3`-end':set(invalid3Clones)},
#                           self.outputDir + self.name +
#                      '_outframe_invalid_primers.png')
#                 del invalid5Clones, invalid3Clones
#             except Exception as e:
#                 raise e
#             del valid3End
#         del OutOfFrame
#         # choose only In-frame RNA sequences
#         self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['v-jframe'] == 'In-frame']
#         # Stop Codon
#         stopcodonInFrameDist = Counter(self.cloneAnnot['stopcodon'].tolist())
#         plotDist(stopcodonInFrameDist, self.name, self.outputDir + self.name +
#                  '_inframe_stopcodon_dist.png', title='In-frame Stop Codons',
#                  proportion=False, rotateLabels=False)
#         print(stopcodonInFrameDist)
#         # stop codon family distribution
#         stopcodFamily = Counter(self.cloneAnnot[self.cloneAnnot['stopcodon'] == 'Yes']['vgene'].tolist())
#         stopcodFamily = compressCountsFamilyLevel(stopcodFamily)
#         plotDist(stopcodFamily, self.name, self.outputDir + self.name +
#                  '_inframe_stopcodfam_dist.png',
#                   title='IGV Abundance of In-frame Unproductive Sequences',
#                  proportion=True)
#         del stopcodonInFrameDist, stopcodFamily
# #         print(stopcodFamily)
#         # choose only productive RNA sequences
#         self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['stopcodon'] == 'No']
#         productiveFamilyDist = compressCountsFamilyLevel(Counter(self.cloneAnnot['vgene'].tolist()))
#         plotDist(productiveFamilyDist, self.name, self.outputDir + self.name +
#                  '_productive_igv_dist.png',
#                   title='IGV Abundance of Productive Sequences',
#                  proportion=True)
#         del productiveFamilyDist
#         if self.end5:
#             valid5End = Counter(self.cloneAnnot['5end'].tolist())
#             plotDist(valid5End, self.name, self.outputDir + self.name +
#                  '_productive_5end_integrity_dist.png', title='Integrity of 5`-end Primer Sequence(Productive)',
#                  proportion=True, rotateLabels=False)
#             invalid5Clones = self.cloneAnnot.index[self.cloneAnnot['5end'] == 'Indelled'].tolist()
#             print("Example of invalid 5`-end:", invalid5Clones[1:10])
#         if self.end3:
#             valid3End = Counter(self.cloneAnnot['3end'].tolist())
#             plotDist(valid3End, self.name, self.outputDir + self.name +
#                  '_productive_3end_integrity_dist.png', title='Integrity of 3`-end Primer Sequence(Productive)',
#                  proportion=True, rotateLabels=False)
#             invalid3Clones = self.cloneAnnot.index[self.cloneAnnot['3end'] == 'Indelled'].tolist()
#             print("Example of invalid 3`-end:", invalid3Clones[1:10])
#             try:
#                 plotVenn({'5`-end':set(invalid5Clones), '3`-end':set(invalid3Clones)},
#                           self.outputDir + self.name +
#                      '_productive_invalid_primers.png')
#             except Exception as e:
#                 raise e
#         gc.collect()
