'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import os
import gzip
import multiprocessing

from collections import Counter, defaultdict

from abseqPy.IgRepAuxiliary.seqUtils import createAlphabet, generateMotif
from abseqPy.IgRepertoire.igRepUtils import writeClonoTypesToFile, createIfNot
from abseqPy.IgRepReporting.igRepPlots import plotSeqLenDist, \
    generateCumulativeLogo, plotSeqDuplication, plotSeqRarefaction, \
    plotSeqRecaptureNew
from abseqPy.logger import LEVEL, printto
from abseqPy.utilities import hasLargeMem, requires


def generateDiversityReport(spectraTypes, clonoTypes, name, outDir, topClonotypes, threads=2, segregate=False,
                            stream=None):
    generateSpectraTypePlots(spectraTypes,  name, outDir, stream=stream)

    flattened = flattenClonoTypeCountsDict(clonoTypes, stream=stream) if segregate else clonoTypes

    writeClonoTypesToFiles(flattened, name, outDir, topClonotypes, stream=stream)
    estimateDiversity(clonoTypes, flattened, name, outDir, threads=threads, segregate=segregate, stream=stream)
#     generateCDRandFRLogos()


def writeClonoTypesToFiles(clonoTypes, name, outDir, topClonotypes=100, stream=None):
    printto(stream, "Clonotype files are being written out ... ")
    cloneFolder = os.path.join(outDir, "clonotypes")

    if not os.path.exists(cloneFolder):
        os.makedirs(cloneFolder)

    for k in clonoTypes.keys():
        # check if the required topClonotypes went overboard, if so, cap to the max length
        if topClonotypes != float('inf') and len(clonoTypes[k]) < topClonotypes:
            stringTopClonotypes = str(len(clonoTypes[k]))
        else:
            stringTopClonotypes = 'all' if topClonotypes == float('inf') else str(topClonotypes)

        # descending order
        filename = os.path.join(cloneFolder, name + ("_{}_clonotypes_{}_over.csv".format(k, stringTopClonotypes)))
        writeClonoTypesToFile(clonoTypes[k], filename, topClonotypes, overRepresented=True)

        # ascending order
        filename = os.path.join(cloneFolder, name + ("_{}_clonotypes_{}_under.csv".format(k, stringTopClonotypes)))
        writeClonoTypesToFile(clonoTypes[k], filename, topClonotypes, overRepresented=False)


def generateSpectraTypePlots(spectraTypes, name, outDir, stream=None):
    specFolder = os.path.join(outDir, "spectratypes")

    if not os.path.exists(specFolder):
        os.makedirs(specFolder)

    for k in spectraTypes.keys():
        filename = os.path.join(specFolder, name + ('_{}_spectratype.csv'.format(k)))
        plotSeqLenDist(spectraTypes[k], name, filename, dna=False,
                       seqName=k.upper(), normed=True, maxbins=40, stream=stream)
        if k == 'cdr3':
            filename = os.path.join(specFolder, name + ('_{}_spectratype_no_outliers.csv'.format(k)))
            plotSeqLenDist(spectraTypes[k], name, filename, dna=False,
                           seqName=k.upper(), normed=True, maxbins=40,
                           removeOutliers=True, stream=stream)


def estimateDiversity(clonoTypes, flatClonoTypes, name, outDir, threads=2, segregate=False, stream=None):
    # create Germline gene level composition logos
    compositionLogos(name, clonoTypes, flatClonoTypes, outDir, threads=threads, detailed=segregate, stream=stream)
    generateSeqMotifs(flatClonoTypes, name, outDir, threads=threads, stream=stream)
    generateRarefactionPlots(flatClonoTypes, name, outDir, threads=threads, stream=stream)
    printto(stream, "The diversity of the library is being estimated ... ")


def generateRarefactionPlots(clonoTypes, name, outDir, threads=2, stream=None):
    regions = clonoTypes.keys()
    regions.sort()
    printto(stream, "Rarefaction files are being generated .... ")
    # select CDR regions only  
    cdrWeights = []
    cdrSeqs = []
    cdrRegions = []
    for region in regions:
        if not region.startswith("cdr"):
            continue
        cdrRegions.append(region.upper())
        cdrSeqs.append(clonoTypes[region].keys())
        cdrWeights.append(map(lambda x: clonoTypes[region][x], cdrSeqs[-1]))
    filename = os.path.join(outDir, name + "_cdr_duplication.csv")
    printto(stream, "\tThe duplication levels is being generated for CDRs .... ")
    plotSeqDuplication(cdrWeights,
                       cdrRegions,
                       filename,
                       'Duplication of CDR Sequences', stream=stream)
    printto(stream, "\tThe rarefaction is being generated for CDRs .... ")
    filename = os.path.join(outDir, name + "_cdr_rarefaction.csv")
    plotSeqRarefaction(cdrSeqs,
                       cdrRegions,
                       filename,
                       cdrWeights,
                       'Rarefaction of CDR Sequences', threads=threads, stream=stream)
    printto(stream, " \tThe percent recapture is being generated for CDRs .... ")
    filename = os.path.join(outDir, name + "_cdr_recapture.csv")
    plotSeqRecaptureNew(cdrSeqs,
                        cdrRegions,
                        filename,
                        'Percent Recapture of CDR Sequences', threads=threads, stream=stream)
    # select FR regions only
    frWeights = []
    frSeqs = []
    frRegions = []
    for region in regions:
        if not region.startswith("fr"):
            continue
        frRegions.append(region.upper())
        frSeqs.append(clonoTypes[region].keys())
        frWeights.append(map(lambda x: clonoTypes[region][x], frSeqs[-1]))
    filename = os.path.join(outDir, name + "_fr_duplication.csv")
    printto(stream, "\tThe duplication levels is being generated for FRs .... ")
    plotSeqDuplication(frWeights,
                       frRegions,
                       filename,
                       'Duplication of FR Sequences', stream=stream)
    printto(stream, "\tThe rarefaction is being generated for FRs .... ")
    filename = os.path.join(outDir,  name + "_fr_rarefaction.csv")
    plotSeqRarefaction(frSeqs,
                       frRegions,
                       filename,
                       frWeights,
                       'Rarefaction of FR Sequences', threads=threads, stream=stream)
    printto(stream, "\tThe percent recapture is being generated for FRs .... ")
    filename = os.path.join(outDir, name + "_fr_recapture.csv")
    plotSeqRecaptureNew(frSeqs,
                        frRegions,
                        filename,
                        'Percent Recapture of FR Sequences', threads=threads, stream=stream)
    # select CDR and V domain 
    cdrWeights = []
    cdrSeqs = []
    cdrRegions = []
    for region in regions:
        if region.startswith("fr"):
            continue
        cdrRegions.append(region.upper())
        cdrSeqs.append(clonoTypes[region].keys())
        cdrWeights.append(map(lambda x: clonoTypes[region][x], cdrSeqs[-1]))
    filename = os.path.join(outDir, name + "_cdr_v_duplication.csv")
    printto(stream, "\tThe duplication levels is being generated for CDRs and V domains .... ")
    plotSeqDuplication(cdrWeights,
                       cdrRegions,
                       filename,
                       'Duplication of CDRs and V Domains', stream=stream)
    printto(stream, "\tThe rarefaction is being generated for CDRs and V domains .... ")
    filename = os.path.join(outDir, name + "_cdr_v_rarefaction.csv")
    plotSeqRarefaction(cdrSeqs,
                       cdrRegions,
                       filename,
                       cdrWeights,
                       'Rarefaction of CDRs and V Domains', threads=threads, stream=stream)
    printto(stream, "\tThe percent recapture is being generated for CDRs and V domains .... ")
    filename = os.path.join(outDir, name + "_cdr_v_recapture.csv")
    plotSeqRecaptureNew(cdrSeqs,
                        cdrRegions,
                        filename,
                        'Percent Recapture of CDRs and V Domains', threads=threads, stream=stream)


def compositionLogos(name, clonoTypes, flatClonoTypes, outDir, threads=2, detailed=False, stream=None):
    """

    :param name: string
                sample name

    :param clonoTypes: dict
    dict with key for each V germline, each having keys of FR / CDR region,
    which in turn, each having a value of Counter() where the AA sequences are tallied
    For example:
    {
        'IGHV3-3': { 'FR1': Counter({"FGWSG": 32, ...}),  'CDR1': Counter(...) },
        'IGHV2-1': { ... }
    }

    :param flatClonoTypes: dict
                    dict with keys of FR / CDR region, each having a value of Counter() where the
                    AA sequences are tallied
                    For example:
                    {
                        'FR1': Counter({"FGWSG": 32, ...}),
                        'CDR1': Counter(...)
                    }

    :param outDir: string

    :param threads: int

    :param detailed: bool
                    segregate composition logo plots based on IGV gene, FR and CDR (all genes combined) composition
                    logos will still be plotted. (If set to false, only FR and CDR composition logos)

    :param stream: stream object
                    output stream
    :return:
    """

    logosFolder = os.path.join(outDir, 'composition_logos')
    createIfNot(logosFolder)
    printto(stream, "Generating composition logos ...")
    if detailed:
        argBuffer = []
        for vgerm in clonoTypes:
            regions = clonoTypes[vgerm].keys()
            regions.sort()
            for region in regions:
                if region == 'v':
                    continue
                clonoType = clonoTypes[vgerm][region]
                seqs = clonoType.keys()
                weights = clonoType.values()

                regionDirectory = os.path.join(logosFolder, region.upper())
                createIfNot(regionDirectory)
                filename = os.path.join(regionDirectory, name + "_{}_cumulative_logo.csv"
                                        .format(vgerm.replace(os.path.sep, '_')))
                if hasLargeMem():
                    printto(stream, "\tbuffering {} for {}".format(region, vgerm))
                    argBuffer.append((seqs, weights, region, filename))
                else:
                    printto(stream, "\tgenerating {} for {}".format(region, vgerm))
                    generateCumulativeLogo(seqs, weights, region, filename, stream=stream)

        if len(argBuffer):
            printto(stream, "Asynchronously generating composition logos from buffer ...")
            pool = multiprocessing.Pool(processes=threads)
            # Generate cumulative sequence logos using Toby's approach
            res = [pool.apply_async(generateCumulativeLogo, args=arg) for arg in argBuffer]
            [p.get() for p in res]  # join processes
            pool.close()
            pool.join()
        printto(stream, "Completed composition logos for IGV families")

    # composition logo for a region(CDR,FR) as a combination of all IGV - i.e. not segregated
    regions = flatClonoTypes.keys()     # combined AA counts(V family) of each region into one sum
    regions.sort()
    for region in regions:
        if region == 'v':
            continue
        clonoType = flatClonoTypes[region]
        seqs = clonoType.keys()
        weights = clonoType.values()

        regionDirectory = os.path.join(logosFolder, region.upper())
        createIfNot(regionDirectory)
        filename = os.path.join(regionDirectory, name + "_cumulative_logo.csv")
        generateCumulativeLogo(seqs, weights, region, filename, stream=stream)


@requires('weblogolib')
def generateSeqMotifs(flatClonoTypes, name, outDir, threads=2, stream=None):
    """
    Create motif plots for FR and CDR regions
    :param flatClonoTypes: dict
                    dict with keys of FR / CDR region, each having a value of Counter() where the
                    AA sequences are tallied
                    For example:
                    {
                        'FR1': Counter({"FGWSG": 32, ...}),
                        'CDR1': Counter(...)
                    }

    :param name: string
                    name of sample

    :param outDir: string
                    output directory

    :param threads: int
                    number of threads to use

    :param stream: stream object
                    output stream
    :return: None
    """

    motifsFolder = os.path.join(outDir, 'motifs')
    createIfNot(motifsFolder)

    printto(stream, "Generating motifs ...")

    # create motif logos
    regions = flatClonoTypes.keys()
    regions.sort()
    argBuffer = []
    for region in regions:
        if region == 'v':
            continue

        clonoType = flatClonoTypes[region]
        seqs = clonoType.keys()
        weights = clonoType.values()

        # Generate sequence motif logos using weblogo
        # generate logos without alignment
        filename = os.path.join(motifsFolder, name + ("_{}_motif_logo.png".format(region)))
        alphabet = createAlphabet(align=False, protein=True, extendAlphabet=True)
        # if hasLargeMem():
        #     printto(stream, "\tbuffering data for {} motif".format(region))
        #     argBuffer.append((seqs, region, alphabet, filename, False, False, True, weights, outDir, threads))
        # else:
        printto(stream, "\tgenerating {} motif".format(region))
        generateMotif(seqs, region, alphabet, filename,  align=False,
                      protein=True, weights=weights, outDir=outDir, threads=threads, stream=stream)

        # generate  logos after alignment
        filename = os.path.join(motifsFolder, name + ("_{}_motif_aligned_logo.png".format(region)))
        alphabet = createAlphabet(align=True, protein=True, extendAlphabet=True)
        # if hasLargeMem():
        #     argBuffer.append((seqs, region, alphabet, filename, True, False, True, weights, outDir, 1))
        # else:
        #     generateMotif(seqs, region, alphabet, filename,  align=True,
        #                   protein=True, weights=weights, outDir=outDir, threads=threads, stream=stream)
        generateMotif(seqs, region, alphabet, filename,  align=True,
                      protein=True, weights=weights, outDir=outDir, threads=threads, stream=stream)

    if len(argBuffer):
        printto(stream, "Asynchronously generating motifs from buffer ...")
        pool = multiprocessing.Pool(processes=threads)
        res = [pool.apply_async(generateMotif, args=arg) for arg in argBuffer]
        [p.get() for p in res]
        pool.close()
        pool.join()

    printto(stream, "CDR/FR Motif analysis complete")


def writeClonotypeDiversityRegionAnalysis(clonoTypes, sampleName, outDir, stream=None):
    """
    For a given set of similar CDR3 clonotypes, it may be classified as a different clonotype if the entire V region
    is considered. This writes the unique counts of other region aside form CDR3s to see if the clonotype will differ
    if the entire V region is considered. Consequently, it's possible to learn which region is (mostly)
    the one responsible of changing the clonotype if it was included.
    :param clonoTypes: DataFrame of clonotypes per read. Requires the CDRs and FRs columns
    :param sampleName: Sample name for output file
    :param outDir: Out directory for output file
    :param stream: debug stream
    :return: None. Produces an output gzipped csv file
    """
    fname = os.path.join(outDir, sampleName + "_clonotype_diversity_region_analysis.csv.gz")
    if os.path.exists(fname):
        printto(stream, "\t File found {}".format(fname), LEVEL.WARN)
        return

    # regions of analysis
    cols = ["cdr1", "cdr2", "fr1", "fr2", "fr3", "fr4"]

    def regionCounts(selectedRows):
        """ returns a list of numbers that corresponds to the frequency of *UNIQUE* "CDR1", "CDR2", .. "FR4"
        (in the order of cols as defined above)
        :param selectedRows: this "DataFrame" of rows should have the same CDR3 region
        :return: a list of numbers, each representing the number of unique region in the order of
        COLS as defined above
        """
        return [str(len(set(selectedRows[region]))) for region in cols]

    # obtain all CDR3s
    cdr3s = set(clonoTypes['cdr3'])

    with gzip.open(fname, "wb") as fp:
        writeBuffer = ""
        # write csv header
        writeBuffer += "cdr3,count," + ','.join(cols) + "\n"
        # for each unique CDR3, find all rows(reads) that have the same CDR3
        for cdr3 in cdr3s:
            rows = clonoTypes[clonoTypes['cdr3'] == cdr3]
            writeBuffer += cdr3 + "," + str(len(rows)) + "," + ','.join(regionCounts(rows)) + '\n'
            if len(writeBuffer) > 4e9:
                fp.write(writeBuffer)
                writeBuffer = ""
        fp.write(writeBuffer)


def flattenClonoTypeCountsDict(clonoTypes, stream=None):
    """
    reduces something of this structure:
            'IGHV1-3': {
                'FR1': { 'FWGCGC': 12, 'EVILK': 1, ... }
                'CDR1': { 'FWGCGC': 12, 'EVILK': 1, ... }
            },
            'IGHV2-3': {
                'FR1' : { 'FWGCGC': 12, 'EVILK': 1, ... }
                'CDR1': { 'FWGCGC': 12, 'EVILK': 1, ... }
            }, ...
    to this:
            {
                'FR1': { 'FWGCGC': 24, 'EVILK': 2, ... }
                'CDR1': { 'FWGCGC': 24, 'EVILK': 2, ... }
            }

    :param clonoTypes: dict
                input nested dictionary

    :return: dict
            flattened dictionary
    """
    printto(stream, "Compressing clonotype table ... discarding IGV information ...")
    flattened = defaultdict(Counter)
    for geneName in clonoTypes:
        for region, counts in clonoTypes[geneName].items():
            flattened[region] += Counter(counts)
    printto(stream, "Finish compressing clonotype table")
    return flattened

