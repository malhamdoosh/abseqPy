'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''
import gc
import os
import logging
import sys

from collections import Counter
from os.path import exists
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas.io.parsers import read_csv
from pandas.io.pytables import read_hdf
from numpy import Inf, random, isnan, logical_not

from abseq.IgRepAuxiliary.upstreamAuxiliary import plotUpstreamLenDist, extractUpstreamSeqs, writeCountsCategoriesToFile
from abseq.IgRepAuxiliary.primerAuxiliary import addPrimerData, generatePrimerPlots
from abseq.config import FASTQC
from abseq.IgRepertoire.igRepUtils import compressCountsGeneLevel, gunzip, fastq2fasta, mergeReads, \
    writeListToFile, writeParams, compressSeqGeneLevel, compressSeqFamilyLevel, loadIGVSeqsFromFasta
from abseq.logger import printto, setupLogger, LEVEL
from abseq.IgRepAuxiliary.productivityAuxiliary import refineClonesAnnotation
from abseq.IgRepReporting.igRepPlots import plotSeqLenDist, plotSeqLenDistClasses, plotVenn, plotDist
from abseq.IgRepAuxiliary.annotateAuxiliary import annotateIGSeqRead
from abseq.IgRepReporting.abundanceReport import writeAbundanceToFiles
from abseq.IgRepReporting.productivityReport import generateProductivityReport
from abseq.IgRepReporting.diversityReport import generateDiversityReport, writeClonotypeDiversityRegionAnalysis
from abseq.IgRepAuxiliary.diversityAuxiliary import annotateSpectratypes, \
    annotateClonotypes
from abseq.IgRepAuxiliary.restrictionAuxiliary import findHits, \
    findHitsRegion, scanRestrictionSitesSimple, loadRestrictionSites
from abseq.IgRepReporting.restrictionReport import generateOverlapFigures


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
    def __init__(self, args):
        self.args = args
        self.task = args.task
        self.reportInterim = args.report_interim
        self.outputDir = args.outdir
        if not os.path.exists(self.outputDir):
            os.makedirs(self.outputDir)
        self.threads = args.threads
        self.primer = args.primer
        self.db = args.database
        self.bitScore = args.bitscore
        self.clonelimit = args.clonelimit
        self.alignLen = args.alignlen
        self.sStart = args.sstart
        self.qStart = args.qstart
        self.seqType = args.seqtype
        self.format = args.fmt
        self.chain = args.chain
        self.name = args.name
        self.fr4cut = args.fr4cut

        if args.task in ['secretion', '5utr']:
            self.upstream = args.upstream

        if args.task in ['rsa', 'rsasimple']:
            self.sitesFile = args.sites

        self.actualQstart = args.actualqstart

        self.trim5End = args.trim5
        self.trim3End = args.trim3

        self.end5 = args.primer5end
        self.end3 = args.primer3end
        self.end5offset = args.primer5endoffset

        self.readFile1 = args.f1
        self.readFile2 = args.f2
        self.merger = args.merger
        self.merge = 'no' if self.merger is None else 'yes'

        self.seqsPerFile = int(10.0 ** 5 / 2)
        self.cloneAnnot = None
        self.cloneSeqs = None
        self.readFile = None
        # True of any of the following directories are already created. We need to distinguish this
        # from the beginning because AbSeq also re-reads HDF within the same analysis to prevent
        # pickling self.cloneAnnot, self.cloneSeqs into multiprocessing.Queue
        self.warnOldDir = any(map(lambda x: exists(self.outputDir + x),
                                  ["abundance/", "productivity/", "diversity/", "restriction_sites/",
                                   "primer_specificity/", 'upstream/']))

        setupLogger(self.name, self.task, args.log)

    def runFastqc(self):
        logger = logging.getLogger(self.name)

        if self.format == 'fasta':
            printto(logger, "Fasta file extension detected, will not perform fastqc", LEVEL.WARN)
            return

        outDir = os.path.join(self.outputDir, "fastqc/")

        if not os.path.isdir(outDir):
            os.makedirs(outDir)

        filename = outDir + self.readFile1.split(os.path.sep)[-1]\
            .replace(".fastq", "")\
            .replace(".gz", "") + "_fastqc.html"

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

    def mergePairedReads(self, outDir=None):
        logger = logging.getLogger(self.name)
        if self.merge != 'yes':
            self.readFile = self.readFile1
        else:
            mergedFastq = mergeReads(self.readFile1, self.readFile2,
                                     self.threads, self.merger, self.outputDir, stream=logger)
            self.readFile = mergedFastq
        if outDir is not None:
            # generate plot of clone sequence length distribution
            outputFile = outDir + self.name + '_all_clones_len_dist.png'
            plotSeqLenDist(self.readFile, self.name, outputFile, self.format,
                           maxbins=40, histtype='bar', removeOutliers=False,
                           normed=True, stream=logger)
            # generate plot of clone sequence length distribution with outliers removed
            outputFile = outDir + self.name + '_all_clones_len_dist_no_outliers.png'
            plotSeqLenDist(self.readFile, self.name, outputFile, self.format,
                           maxbins=40, histtype='bar', removeOutliers=True,
                           normed=True, stream=logger)

    def annotateClones(self, outDirFilter=None, all=False):
        logger = logging.getLogger(self.name)
        outDir = os.path.join(self.outputDir, "annot/")

        if not os.path.isdir(outDir):
            os.makedirs(outDir)

        cloneAnnotFile = outDir + self.name + "_clones_annot.h5"

        if self.readFile is None:
            self.mergePairedReads(outDir)

        if exists(cloneAnnotFile):
            if self.task == "annotate":
                printto(logger, "\tClones annotation file found and no further work needed ... " +
                        cloneAnnotFile.split('/')[-1])
            else:
                printto(logger, "\tClones annotation file found and being loaded ... " +
                        cloneAnnotFile.split('/')[-1])
                self.cloneAnnot = read_hdf(cloneAnnotFile, "cloneAnnot")
        else:
            if not exists(self.readFile):
                raise Exception(self.readFile + " does not exist!")

            # Convert FASTQ file into FASTA format
            if self.format == 'fastq':
                readFasta = fastq2fasta(self.readFile, self.outputDir, stream=logger)
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
                                                               self.seqType, outdir=outDir, stream=logger)
            sys.stdout.flush()
            gc.collect()

            if (len(filteredIDs) > 0):
                writeListToFile(filteredIDs, outDir + self.name + "_unmapped_clones.txt")
            # export the CDR/FR annotation to a file
            printto(logger, "\tClones annotation file is being written to " +
                    cloneAnnotFile.split("/")[-1])
            #             self.cloneAnnot.to_csv(cloneAnnotFile, sep='\t', header=True, index=True)
            self.cloneAnnot.to_hdf(cloneAnnotFile, "cloneAnnot", mode='w')
            paramFile = writeParams(self.args, outDir)
            printto(logger, "The analysis parameters have been written to " + paramFile)

        printto(logger, "Number of clones that are annotated is {0:,}".format(
                int(self.cloneAnnot.shape[0])), LEVEL.INFO)
        if outDirFilter or all:
            if outDirFilter is None:
                outDirFilter = outDir
            # Filter clones based on bitscore, alignLen, qStart, and sStart
            printto(logger, "Clones are being filtered based on the following criteria: ", LEVEL.INFO)
            printto(logger, "\tBit score: " + repr(self.bitScore), LEVEL.INFO)
            printto(logger, "\tAlignment length: " + repr(self.alignLen), LEVEL.INFO)
            printto(logger, "\tSubject V gene start: " + repr(self.sStart), LEVEL.INFO)
            printto(logger, "\tQuery V gene start: " + repr(self.qStart), LEVEL.INFO)
            selectedRows = (
                    (self.cloneAnnot['bitscore'] >= self.bitScore[0]) &  # check bit-Score
                    (self.cloneAnnot['bitscore'] <= self.bitScore[1]) &
                    (self.cloneAnnot['alignlen'] >= self.alignLen[0]) &  # check alignment length
                    (self.cloneAnnot['alignlen'] <= self.alignLen[1]) &
                    (self.cloneAnnot['vstart'] >= self.sStart[0]) &  # check subject (V gene) start position
                    (self.cloneAnnot['vstart'] <= self.sStart[1]) &
                    (self.cloneAnnot['vqstart'] >= self.qStart[0]) &  # check query (V gene) start position
                    (self.cloneAnnot['vqstart'] <= self.qStart[1])
            )
            filteredIDs = self.cloneAnnot[logical_not(selectedRows)]
            if len(filteredIDs) > 0:
                filteredIDs = filteredIDs[['vgene', 'vstart', 'vqstart', 'bitscore', 'alignlen']]
                filteredIDs.to_csv(outDirFilter + self.name + "_filtered_out_clones.txt",
                                   sep="\t", header=True, index=True)
            retained = len(selectedRows) - len(filteredIDs)
            printto(logger, 'Percentage of retained clones is {0:,.2f}% ({1:,}/{2:,})'.format(
                    retained * 100.0 / self.cloneAnnot.shape[0],
                    retained,
                    int(self.cloneAnnot.shape[0]), LEVEL.INFO))
            self.cloneAnnot = self.cloneAnnot[selectedRows]

    def analyzeAbundance(self, all=False):
        # Estimate the IGV family abundance for each library
        logger = logging.getLogger(self.name)
        outDir = os.path.join(self.outputDir, "abundance/")

        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: remove the 'abundance' directory if you changed the filtering criteria.",
                    LEVEL.WARN)
        if not all or self.cloneAnnot is None:
            self.annotateClones(outDir)

        writeAbundanceToFiles(self.cloneAnnot, self.name, outDir, self.chain, stream=logger)
        gc.collect()
        paramFile = writeParams(self.args, outDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzeProductivity(self, generateReport=True, all=False, inplace=True):
        """
        Analyzes productivity
        :param generateReport:
        :param all:
        :param inplace: NOTE: if this is set to true, self.cloneAnnot and self.cloneSeqs will only contain
        productive sequences after this method finishes. Set to false to retain all sequences
        :return:
        """
        logger = logging.getLogger(self.name)
        outDir = os.path.join(self.outputDir, "productivity/")

        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: remove the 'productivity' directory if you changed the filtering criteria.",
                    LEVEL.WARN)

        refinedCloneAnnotFile = outDir + self.name
        refinedCloneAnnotFile += "_refined_clones_annot.h5"

        cloneSeqFile = outDir + self.name
        cloneSeqFile += "_clones_seq.h5"

        if not exists(refinedCloneAnnotFile):
            if not all or self.cloneAnnot is None:
                self.annotateClones(outDir)
            #             if (self.trimmed):
            #                 self.trim3End = 0
            #                 self.trim5End = 0
            #             elif self.trim3End > 0 or self.trim5End > 0:
            #                 print("WARNING: if trimming was applied in the 'annotate' step, you may not need trimming")
            # print(sys.getsizeof(self.cloneAnnot) / (1024.**3)) # in GB
            (self.cloneAnnot, self.cloneSeqs) = refineClonesAnnotation(outDir, self.name,
                                                                       self.cloneAnnot, self.readFile,
                                                                       self.format, self.actualQstart,
                                                                       self.chain, self.fr4cut,
                                                                       self.trim5End, self.trim3End,
                                                                       self.seqsPerFile, self.threads, self.db,
                                                                       stream=logger)
            gc.collect()
            # if generateReport:
            # export the CDR/FR annotation to a file                
            printto(logger, "The refined clone annotation file is being written to " +
                    refinedCloneAnnotFile.split("/")[-1])
            self.cloneAnnot.to_hdf(refinedCloneAnnotFile, "refinedCloneAnnot", mode='w',
                                   complib='blosc')
            sys.stdout.flush()
            printto(logger, "The clone protein sequences are being written to " +
                    cloneSeqFile.split("/")[-1])
            self.cloneSeqs.to_hdf(cloneSeqFile, "cloneSequences", mode='w',
                                  complib='blosc')
            sys.stdout.flush()
            paramFile = writeParams(self.args, outDir)
            printto(logger, "The analysis parameters have been written to " + paramFile)
        else:
            printto(logger, "The refined clone annotation files were found and being loaded ... " +
                    refinedCloneAnnotFile.split('/')[-1])
            self.cloneAnnot = read_hdf(refinedCloneAnnotFile, "refinedCloneAnnot")
            printto(logger, "\tClone annotation was loaded successfully")
            self.cloneSeqs = read_hdf(cloneSeqFile, "cloneSequences")
            printto(logger, "\tClone sequences were loaded successfully")
        if generateReport:
            # display statistics
            printto(logger, "Productivity report is being generated ... ")
            generateProductivityReport(self.cloneAnnot, self.cloneSeqs, self.name, self.chain, outDir, stream=logger)
            # TODO: analyze productive clones only
        #             self.analyzeIgProtein()
        #             sys.stdout.flush()
        # Diversity analysis can be applied on productive clones only     
        before = int(self.cloneAnnot.shape[0])
        inFrame = self.cloneAnnot[self.cloneAnnot['v-jframe'] == 'In-frame']
        if inplace:
            cloneAnnot = self.cloneAnnot = inFrame[inFrame['stopcodon'] == 'No']
            self.cloneSeqs = self.cloneSeqs.loc[self.cloneAnnot.index]
        else:
            cloneAnnot = inFrame[inFrame['stopcodon'] == 'No']
        printto(logger, "Percentage of productive clones {0:,.2f}% ({1:,}/{2:,})".format(
            cloneAnnot.shape[0] * 100.0 / before,
            int(cloneAnnot.shape[0]),
            int(before)
            ), LEVEL.INFO)

    '''
    
    '''

    def analyzeDiversity(self, all=False):
        logger = logging.getLogger(self.name)
        outDir = os.path.join(self.outputDir,  "diversity/")

        if not all or self.cloneAnnot is None or self.cloneSeqs is None:
            self.analyzeProductivity(self.reportInterim, all)

        if len(self.cloneAnnot) == 0:
            printto(logger, "WARNING: There are no productive sequences found (post-refinement) in {},"
                            " skipping diversity analysis.".format(self.name), LEVEL.WARN)
            return

        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: remove the 'diversity' directory if you changed the filtering criteria.",
                    LEVEL.WARN)

        gc.collect()

        # Identify spectratypes 
        printto(logger, "Spectratypes are being calculated ... ")
        spectraTypes = annotateSpectratypes(self.cloneAnnot)

        # Identify clonotypes 
        printto(logger, "Clonotypes are being generated ... ")
        clonoTypes = annotateClonotypes(self.cloneSeqs)

        generateDiversityReport(spectraTypes, clonoTypes, self.name, outDir,
                                self.clonelimit, stream=logger)

        # todo: remove this for now - it's unoptimized and extremely slow
        # writeClonotypeDiversityRegionAnalysis(self.cloneSeqs, self.name, outDir, stream=logger)

        paramFile = writeParams(self.args, outDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzeRestrictionSitesSimple(self):
        # TODO: parallelize this function to run faster
        logger = logging.getLogger(self.name)

        outDir = self.outputDir + "restriction_sites/"

        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        elif self.warnOldDir:
            print("WARNING: remove the 'restriction_sites' directory if you changed the filtering criteria.")

        siteHitsFile = outDir + self.name
        siteHitsFile += "_%s_rsasimple.csv" % (self.sitesFile.split('/')[-1].split('.')[0])
        overlap2File = siteHitsFile.replace('.csv', '_overlap_order2.csv')

        if (exists(siteHitsFile)):
            print("Restriction sites were already scanned at ... " +
                  siteHitsFile.split('/')[-1])
            rsaResults = read_csv(siteHitsFile, header=0)
            if exists(overlap2File):
                overlapResults = {}
                overlapResults['order2'] = read_csv(overlap2File, header=0, index_col=0)
            else:
                overlapResults = None
        else:
            self.annotateClones(outDir)
            sys.stdout.flush()
            (rsaResults, overlapResults) = scanRestrictionSitesSimple(self.name,
                                                                      self.readFile, self.format,
                                                                      self.cloneAnnot, self.sitesFile,
                                                                      self.threads)
            rsaResults.to_csv(siteHitsFile,
                              header=True,
                              index=False)
            print("RSA results were written to " + siteHitsFile.split("/")[-1])
            if overlapResults.get("order2", None) is not None:
                overlapResults["order2"].to_csv(overlap2File,
                                                header=True, index=True)
        # # print out the results        
        generateOverlapFigures(overlapResults,
                               rsaResults.loc[rsaResults.shape[0] - 1, "No.Molecules"],
                               self.name, siteHitsFile)
        paramFile = writeParams(self.args, outDir)
        printto(logger, "The analysis parameters have been written to " + paramFile)

    def analyzeRestrictionSites(self, all=False):
        #         sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'
        #         sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        siteHitsFile = self.outputDir + self.name
        siteHitsFile += "_%s.csv" % (self.sitesFile.split('/')[-1].split('.')[0])

        if exists(siteHitsFile):
            print("Restriction sites were already searched at ... " + siteHitsFile.split('/')[-1])
            return

        if not all or self.cloneAnnot is None or self.cloneSeqs is None:
            self.analyzeProductivity(self.reportInterim, all)

        rsites = loadRestrictionSites(self.sitesFile)
        print("Restriction sites are being searched ... ")
        gc.collect()
        self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['v-jframe'] == 'In-frame']
        self.cloneAnnot = self.cloneAnnot[self.cloneAnnot['stopcodon'] == 'No']
        queryIds = self.cloneAnnot.index
        siteHitsCount = {}
        siteHitSeqsCount = {}
        hitRegion = {}
        siteHitSeqsGermline = {}
        seqsCutByAny = 0
        siteHitsSeqsIDs = {}
        siteHitsSeqsIGV = {}
        for site in rsites.keys():
            siteHitsCount[site] = 0
            siteHitSeqsCount[site] = 0
            hitRegion[site] = Counter({'fr1': 0, 'cdr1': 0,
                                       'fr2': 0, 'cdr2': 0,
                                       'fr3': 0, 'cdr3': 0,
                                       'fr4': 0})
            siteHitSeqsGermline[site] = []
            siteHitsSeqsIDs[site] = set()
            siteHitsSeqsIGV[site] = set()
        germline = {'fr1', 'fr2', 'fr3', 'cdr1', 'cdr2'}
        procSeqs = 0
        #         if (MEM_GB > 20):
        #             TODO: remember to make sure SeqIO.parse is parsing a unzipped self.readFile1
        #                   (use safeOpen from IgRepertoire.utils) if not sure
        #             records = SeqIO.to_dict(SeqIO.parse(self.readFile1, self.format))
        #         else:
        # SeqIO.index can only open string file names and they must be uncompressed
        records = SeqIO.index(gunzip(self.readFile1), self.format)
        for id in queryIds:
            record = records[id]
            try:
                qsRec = self.cloneAnnot.loc[record.id].to_dict()
                qstart = qsRec['vqstart'] - qsRec['vstart']  # zero-based
                if (isnan(qsRec['fr4.end'])):
                    end = len(record.seq)
                else:
                    end = int(qsRec['fr4.end'])
                seq = str(record.seq[qstart:end])
                seqRC = str(Seq(seq).reverse_complement())
                cut = False
                for site in siteHitsCount.keys():
                    hits = findHits(seq, rsites[site])
                    strand = "forward"
                    if len(hits) == 0:
                        hits = findHits(seqRC, rsites[site])
                        strand = "reversed"
                    if len(hits) > 0:
                        siteHitsCount[site] += len(hits)
                        siteHitSeqsCount[site] += 1
                        hitsRegion = findHitsRegion(qsRec, hits)
                        if (len(set(hitsRegion).intersection(germline)) > 0
                                and len(siteHitSeqsGermline[site]) < 10000):
                            siteHitSeqsGermline[site].append((strand, str(record.seq)))
                            siteHitsSeqsIGV[site].add(qsRec['vgene'].split('*')[0])
                        hitRegion[site] += Counter(hitsRegion)
                        siteHitsSeqsIDs[site].add(record.id)
                        cut = True
                if cut:
                    seqsCutByAny += 1
                procSeqs += 1
                if procSeqs % self.seqsPerFile == 0:
                    print('%d/%d sequences have been searched ... ' % (procSeqs, len(queryIds)))
            #                 break
            except BaseException as e:
                print(qstart, end, len(record.seq), str(record.seq))
                print(e)
                raise
        records.close()
        print('%d/%d sequences have been searched ... ' % (procSeqs, len(queryIds)))
        # # print out the results
        f = open(siteHitsFile, 'w')
        f.write("Enzyme,Restriction Site,No.Hits,Percentage of Hits (%),"
                "No.Molecules,Percentage of Molecules (%),FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4, V Germlines \n")
        sites = sorted(siteHitSeqsCount, key=siteHitSeqsCount.get)
        for site in sites:
            f.write("%s,%s,%d,%.3f,%d,%.3f,%d,%d,%d,%d,%d,%d,%d,%s\n" % (site,
                                                                         rsites[site],
                                                                         siteHitsCount[site],
                                                                         siteHitsCount[site] * 100.0 / sum(
                                                                             siteHitsCount.values()),
                                                                         siteHitSeqsCount[site],
                                                                         siteHitSeqsCount[site] * 100.0 / len(queryIds),
                                                                         hitRegion[site]['fr1'],
                                                                         hitRegion[site]['cdr1'],
                                                                         hitRegion[site]['fr2'],
                                                                         hitRegion[site]['cdr2'],
                                                                         hitRegion[site]['fr3'],
                                                                         hitRegion[site]['cdr3'],
                                                                         hitRegion[site]['fr4'],
                                                                         '|'.join(siteHitsSeqsIGV[site])))
            # write the first 100 sequences cut in the germline of each restriction enzyme
            seqs = []
            for (strand, seq) in siteHitSeqsGermline[site]:
                seqs.append(SeqRecord(Seq(seq), id='seq' + `len(seqs)` + strand))
            SeqIO.write(seqs, siteHitsFile.replace('.csv', '_germline' + site + '.fasta'), 'fasta')
        f.write("Sequences cut by any of the above enzymes, %d, %.3f\n" % (seqsCutByAny,
                                                                           seqsCutByAny * 100.0 / len(queryIds)))
        f.close()
        # Ven Diagram of overlapping sequences
        plotVenn(siteHitsSeqsIDs, siteHitsFile.replace('.csv', '_venn.png'))
        print("Restriction enzyme results were written to " + siteHitsFile)

    def analyzeSecretionSignal(self, all=False):
        logger = logging.getLogger(self.name)
        outputDir = os.path.join(self.outputDir, 'upstream/')

        if os.path.exists(outputDir):
            os.makedirs(outputDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: Remove 'upstream' directory if you've changed the filtering criteria.",
                    LEVEL.WARN)

        # need self.cloneAnnot dataframe for further analysis
        if self.cloneAnnot is None or not all:
            self.analyzeAbundance(all=all)

        printto(logger, "The diversity of the upstream of IGV genes is being analyzed ... ")

        upstreamFile = outputDir + self.name + "_upigv_{:.0f}_{:.0f}.fasta"\
            .format(self.upstream[0], self.upstream[1])

        if not exists(upstreamFile):
            extractUpstreamSeqs(self.cloneAnnot, self.readFile, self.upstream, upstreamFile, stream=logger)
        else:
            printto(logger, "\tUpstream sequences file was found! ... " + os.path.basename(upstreamFile), LEVEL.WARN)

        upstreamFile = os.path.abspath(upstreamFile)

        expectLength = self.upstream[1] - self.upstream[0] + 1

        # plot the distribution of sequence length
        plotUpstreamLenDist(upstreamFile, expectLength, self.name, stream=logger)

        if expectLength != Inf:
            # classify secretion signals based on length, ATG location, gene and gene family
            # analyze intact secretion signals
            self.analyzeSequences(upstreamFile, self.name, [expectLength, expectLength], True)
            # analyze trimmed secretion signals
            self.analyzeSequences(upstreamFile, self.name, [1, expectLength - 1], True)

    def analyze5UTR(self, all=False):
        logger = logging.getLogger(self.name)
        outputDir = os.path.join(self.outputDir, 'upstream/')

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: Remove 'upstream' directory if you've changed the filtering criteria",
                    LEVEL.WARN)

        # requires self.cloneAnnot dataframe for further analysis
        if self.cloneAnnot is None or not all:
            self.analyzeAbundance(all=all)

        printto(logger, "The diversity of the upstream of IGV genes is being analyzed ... ")

        upstreamFile = outputDir + self.name + "_5utr_{:.0f}_{:.0f}.fasta"\
            .format(self.upstream[0], self.upstream[1])

        if not exists(upstreamFile):
            extractUpstreamSeqs(self.cloneAnnot, self.readFile, self.upstream, upstreamFile, stream=logger)
        else:
            printto(logger, "\tUpstream sequences file was found! ... " + os.path.basename(upstreamFile), LEVEL.WARN)

        upstreamFile = os.path.abspath(upstreamFile)
        expectLength = self.upstream[1] - self.upstream[0] + 1

        plotUpstreamLenDist(upstreamFile, expectLength, self.name, stream=logger)

        # if user provided values to upstream (and it's not Inf)
        if expectLength != Inf:
            # analyze intact secretion signals
            self.analyzeSequences(upstreamFile,
                                  self.name,
                                  [expectLength, expectLength],
                                  startCodon=False,
                                  type='5utr',
                                  clusterMotifs=True)

    def analyzePrimerSpecificity(self, all=False):
        logger = logging.getLogger(self.name)
        outDir = os.path.join(self.outputDir, 'primer_specificity/')

        if not os.path.exists(outDir):
            os.makedirs(outDir)
        elif self.warnOldDir:
            printto(logger, "WARNING: remove the 'primer_specificity' directory if "
                            "you changed the filtering criteria.", LEVEL.WARN)

        # the naming convention obeys previous h5 dataframes:
        # if self.cloneAnnot was refined, primer_annot will be named samplename_primer_annot_refined.h5
        # else if self.cloneAnnot was not refined, primer_annot will be named samplename_primer_annot.h5
        primerAnnotFile = os.path.join(outDir, self.name + "_primer_annot{}.h5"
                                       .format('_refined' if self.fr4cut else ''))

        # if we can't find hdf file, create it, else read it
        if not exists(primerAnnotFile):
            # Load self.cloneAnnot for further analysis.
            # skip checking for existence of dataframes, analyzeProd/Abun will do it for us
            if not all or self.cloneAnnot is None:
                if exists(os.path.join(self.outputDir, 'productivity',
                                       self.name + '_refined_clones_annot.h5')):
                    printto(logger, "Using refined clone annotation for primer specificity analysis")
                    self.analyzeProductivity(all=all, inplace=False)
                else:
                    printto(logger, "Using unrefined clone annotation for primer specificity analysis")
                    self.analyzeAbundance(all)
            # add additional primer related data to the dataframe generated by either abundance/productivity analysis
            # before we begin primer analysis
            self.cloneAnnot = addPrimerData(self.cloneAnnot, self.readFile, self.format, self.fr4cut,
                                            self.trim5End, self.trim3End, self.actualQstart,
                                            self.end5, self.end3, self.end5offset, self.threads, stream=logger)
            # save new "extended dataframe" into primer_specificity directory
            self.cloneAnnot.to_hdf(primerAnnotFile,
                                   "primerCloneAnnot", mode='w', complib='blosc')

        else:
            printto(logger, "The{}primer clone annotation files were found and being loaded ... "
                  .format(' refined ' if self.fr4cut else ' '), LEVEL.WARN)
            self.cloneAnnot = read_hdf(primerAnnotFile, "primerCloneAnnot")
            printto(logger, "\tPrimer clone annotation was loaded successfully")

        # TODO: Fri Feb 23 17:13:09 AEDT 2018
        # TODO: check findBestMatchAlignment of primer specificity best match, see if align.localxx is used correctly!
        generatePrimerPlots(self.cloneAnnot, outDir, self.name, self.end5, self.end3, stream=logger)

    def analyzeSeqLen(self, klass=False):
        logger = logging.getLogger(self.name)
        outdir = os.path.join(self.args.outdir, 'annot/')

        printto(logger, "Sequence {}length distribution is being calculated ... ".format('class' if klass else ''))

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if klass:
            outputFile = outdir + self.args.name + '_length_dist_classes.png'
            plotSeqLenDistClasses(self.args.f1, self.args.name, outputFile, self.args.fmt, stream=logger)
        else:
            outputFile = outdir + self.args.name + '_seq_length_dist.png'
            plotSeqLenDist(self.args.f1, self.args.name, outputFile, self.args.fmt, maxbins=-1, stream=logger)

    def loadValidSequences(self, upstreamFile, sampleName, expectLength, startCodon=True, type='secsig'):
        print("\tSequences between %d and %d are being extracted ... "
              % (expectLength[0], expectLength[1]))
        ighvSignals = {}
        ighvSignalsCounts = Counter()
        ighvSignalsNoATG = []
        noStartCodonCounts = Counter()
        faultyTrans = []
        faultyTransCounts = Counter()
        #         if (MEM_GB > 20):
        #             TODO: remember to make sure SeqIO.parse is parsing a unzipped self.readFile1
        #                   (use safeOpen from IgRepertoire.utils) if not sure
        #             records = SeqIO.to_dict(SeqIO.parse(upstreamFile, 'fasta'))
        #         else:
        # SeqIO.index can only parse string filename (that isn't opened) and unzipped
        records = SeqIO.index(gunzip(upstreamFile), 'fasta')
        for id in records:
            rec = records[id]
            ighv = rec.id.split('|')[1]
            seq = rec.seq
            if (expectLength[0] <= len(rec) <= expectLength[1]):
                if (not startCodon or "ATG" in seq):
                    if (faultyTransCounts.get(ighv, None) is None):
                        faultyTransCounts[ighv] = 0
                    if (type == 'secsig'):
                        seq = seq.translate(to_stop=False)[1:]
                    if ('X' in seq or '*' in seq):
                        # print(rec.id, str(rec.seq), str(seq))
                        faultyTrans.append(rec)
                        faultyTransCounts[ighv] += 1
                    elif ('N' not in rec.seq):
                        # rec.seq = seq
                        if (ighvSignals.get(ighv, None) is None):
                            ighvSignals[ighv] = []
                            ighvSignalsCounts[ighv] = 0
                        ighvSignals[ighv].append(str(rec.seq))  # rec
                        ighvSignalsCounts[ighv] += 1
                    else:
                        print('Ignored: ' + str(rec.seq) + ' ' + str(seq))
                        if (type == 'secsig'):
                            faultyTrans.append(rec)
                            faultyTransCounts[ighv] += 1
                elif startCodon:
                    ighvSignalsNoATG.append(rec)  # seq
                    if (noStartCodonCounts.get(ighv, None) is None):
                        noStartCodonCounts[ighv] = 0
                    noStartCodonCounts[ighv] += 1
        records.close()
        if (sum(ighvSignalsCounts.values()) > 0):
            print("\tThere are %d VALID secretion signals within expected length %s and startCodon=%s " %
                  (sum(ighvSignalsCounts.values()), str(expectLength), startCodon))
            if (type == 'secsig'):
                title = 'Valid Secretion Signals'
            else:
                title = 'Valid 5`-UTRs'
            writeCountsCategoriesToFile(ighvSignalsCounts, sampleName,
                                        self.outputDir + sampleName + '_%s%d%d_valid_' % (
                                        type, expectLength[0], expectLength[1]),
                                        title)
        # # Faulty secretion signals: stop codons or low quality sequencing
        if (len(faultyTrans) > 0):
            # variant level
            faultySeqFile = self.outputDir + sampleName + '_%s%d%d_faulty_trans.fasta' % (
            type, expectLength[0], expectLength[1])
            SeqIO.write(faultyTrans, faultySeqFile, 'fasta')
            writeCountsCategoriesToFile(faultyTransCounts, sampleName,
                                        self.outputDir + sampleName + '_%s%d%d_faulty_' % (
                                        type, expectLength[0], expectLength[1]),
                                        'Faulty Translations')
            print("\tTotal faulty secretion signals is %d (excluded)" % (len(faultyTrans)))
            examples = random.choice(range(len(faultyTrans)), 5)
            for i in examples:
                print(faultyTrans[i].seq, faultyTrans[i].seq.translate())
            faultyTrans = None
        else:
            faultySeqFile = None
            # secretion signals with no start codons
        if (len(ighvSignalsNoATG) > 0):
            noStartCodonFile = self.outputDir + sampleName + '_%s%d%d_no_atg.fasta' % (
            type, expectLength[0], expectLength[1])
            SeqIO.write(ighvSignalsNoATG, noStartCodonFile, 'fasta')
            writeCountsCategoriesToFile(noStartCodonCounts, sampleName,
                                        self.outputDir + sampleName + '_%s%d%d_no_atg_' % (
                                        type, expectLength[0], expectLength[1]),
                                        'Secretion Signals without Start Codon')
            print("\tThere is no ATG codon in %d sequences (excluded). " % (len(ighvSignalsNoATG)))
            examples = random.choice(range(len(ighvSignalsNoATG)), 5)
            for i in examples:
                print(ighvSignalsNoATG[i].seq)
            ighvSignalsNoATG = None
        else:
            noStartCodonFile = None
        gc.collect()
        return (ighvSignals, faultySeqFile, noStartCodonFile)

    def analyzeSequences(self, upstreamFile, sampleName, expectLength, startCodon=True,
                         type='secsig', clusterMotifs=False):
        from abseq.IgRepAuxiliary.SeqUtils import generateMotifs
        lastFile = self.outputDir + sampleName + '_%s%d%d_dna_family' % (type, expectLength[0], expectLength[1])
        lastFile += '_consensus.txt'
        if (exists(lastFile)):
            print("Sequences were already analyzed " + lastFile)
            ighvSignals = {}
            faultySeqFile = self.outputDir + sampleName + '_%s%d%d_faulty_trans.fasta' % (
                type, expectLength[0], expectLength[1])
            noStartCodonFile = self.outputDir + sampleName + '_%s%d%d_no_atg.fasta' % (
                type, expectLength[0], expectLength[1])
        else:
            print("Sequences are being analyzed ... ")
            (ighvSignals, faultySeqFile, noStartCodonFile) = self.loadValidSequences(upstreamFile,
                sampleName, expectLength, startCodon, type)

        # extract DNA motifs for each germline variant
        generateMotifs(ighvSignals, expectLength[0] < expectLength[1],
                       self.outputDir + sampleName +
                       '_%s%d%d_dna_variant' % (type, expectLength[0], expectLength[1]),
                       clusterMotifs=clusterMotifs)
        # extract protein motifs for each each germline variant
        if expectLength[0] == expectLength[1] and type == 'secsig':
            faultySeq = loadIGVSeqsFromFasta(faultySeqFile)
            generateMotifs(faultySeq, True,
                           self.outputDir + sampleName +
                           '_%s%d%d_faulty_variant' % (type, expectLength[0], expectLength[1]),
                           transSeq=False, extendAlphabet=True,
                           clusterMotifs=clusterMotifs)
            noStartCodonSeq = loadIGVSeqsFromFasta(noStartCodonFile)
            generateMotifs(noStartCodonSeq, True,
                           self.outputDir + sampleName +
                           '_%s%d%d_untranslated_variant' % (type, expectLength[0], expectLength[1]),
                           transSeq=False, extendAlphabet=True,
                           clusterMotifs=clusterMotifs)
            generateMotifs(ighvSignals, False,
                           self.outputDir + sampleName + '_%s%d%d_protein_variant' % (
                           type, expectLength[0], expectLength[1]),
                           transSeq=True,
                           clusterMotifs=clusterMotifs)

            # extract motifs for germline genes
        ighvSignals = compressSeqGeneLevel(ighvSignals)
        generateMotifs(ighvSignals, expectLength[0] < expectLength[1],
                       self.outputDir + sampleName + '_%s%d%d_dna_gene' % (type, expectLength[0], expectLength[1]),
                       clusterMotifs=clusterMotifs)
        if expectLength[0] == expectLength[1] and type == 'secsig':
            faultySeq = compressSeqGeneLevel(faultySeq)
            generateMotifs(faultySeq, True,
                           self.outputDir + sampleName +
                           '_%s%d%d_faulty_gene' % (type, expectLength[0], expectLength[1]),
                           transSeq=False, extendAlphabet=True,
                           clusterMotifs=clusterMotifs)
            noStartCodonSeq = compressSeqGeneLevel(noStartCodonSeq)
            generateMotifs(noStartCodonSeq, True,
                           self.outputDir + sampleName +
                           '_%s%d%d_untranslated_gene' % (type, expectLength[0], expectLength[1]),
                           transSeq=False, extendAlphabet=True,
                           clusterMotifs=clusterMotifs)
            generateMotifs(ighvSignals, False,
                           self.outputDir + sampleName + '_%s%d%d_protein_gene' % (
                           type, expectLength[0], expectLength[1]),
                           transSeq=True,
                           clusterMotifs=clusterMotifs)

            # extract motifs for germline families
        ighvSignals = compressSeqFamilyLevel(ighvSignals)
        generateMotifs(ighvSignals, expectLength[0] < expectLength[1],
                       self.outputDir + sampleName + '_%s%d%d_dna_family' % (type, expectLength[0], expectLength[1]),
                       clusterMotifs=clusterMotifs)
        if expectLength[0] == expectLength[1] and type == 'secsig':
            faultySeq = compressSeqFamilyLevel(faultySeq)
            generateMotifs(faultySeq, True,
                           self.outputDir + sampleName +
                           '_%s%d%d_faulty_family' % (type, expectLength[0], expectLength[1]),
                           transSeq=False, extendAlphabet=True,
                           clusterMotifs=clusterMotifs)
            noStartCodonSeq = compressSeqFamilyLevel(noStartCodonSeq)
            generateMotifs(noStartCodonSeq, True,
                           self.outputDir + sampleName +
                           '_%s%d%d_untranslated_family' % (type, expectLength[0], expectLength[1]),
                           transSeq=False, extendAlphabet=True,
                           clusterMotifs=clusterMotifs)
            generateMotifs(ighvSignals, False,
                           self.outputDir + sampleName + '_%s%d%d_protein_family' % (
                           type, expectLength[0], expectLength[1]),
                           transSeq=True,
                           clusterMotifs=clusterMotifs)

    def analyzeIgProtein(self):
        # sampleName = self.readFile1.split('/')[-1].split("_")[0] + '_'
        # sampleName += self.readFile1.split('/')[-1].split("_")[-1].split('.')[0]
        # self.cloneSeqs = read_csv(cloneSeqFile, sep='\t',
        #                                header=0, index_col=0)
        self.readFile1 = self.outputDir + self.name
        self.readFile1 += '_productive_prot.fasta'
        if (not exists(self.readFile1)):
            print("Protein sequences are being prepared ...")
            records = []
            procSeqs = 0
            open(self.readFile1, 'w').close()
            for id in self.cloneSeqs.index:
                seq = ''.join(self.cloneSeqs.loc[id,].tolist()[1:])
                if '*' in seq:
                    seq = seq.replace('*', 'X')
                rec = SeqRecord(Seq(seq), id=id, name="", description="")
                records.append(rec)
                procSeqs += 1
                if procSeqs % self.seqsPerFile == 0:
                    print('\t%d/%d sequences have been processed ...  ' % (procSeqs, len(self.cloneSeqs.index)))
                    sys.stdout.flush()
                    SeqIO.write(records, open(self.readFile1, 'a'), 'fasta')
                    records = []
            SeqIO.write(records, open(self.readFile1, 'a'), 'fasta')
            del records
        else:
            print("File found ... " + self.readFile1.split('/')[-1])
        self.format = 'fasta'
        self.readFile2 = None
        self.seqType = 'protein'
        self.bitScore = [0, Inf]
        self.alignLen = [0, Inf]
        self.sStart = [1, Inf]
        if (exists(self.outputDir + "/abundance/")):
            print("Protein sequences have been already analyzed ... ")
        else:
            self.analyzeAbundance()


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
