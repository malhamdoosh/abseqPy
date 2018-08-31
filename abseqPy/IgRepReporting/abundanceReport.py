'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
'''

import os
from collections import Counter, defaultdict

from abseqPy.IgRepReporting.igRepPlots import plotDist, generateStatsHeatmap, writeCSV
from abseqPy.IgRepertoire.igRepUtils import compressCountsGeneLevel, \
    compressCountsFamilyLevel
from abseqPy.logger import printto, LEVEL


def writeVAbundanceToFiles(stats, sampleName, outDir, stream=None):
    igvDist = Counter(stats["vgene"].tolist())
    if len(igvDist) == 0:
        printto(stream, "WARNING: No IGV hits were detected.", LEVEL.WARN)
        return

    # Write the counts of all IGVs into a text file - variant_level isn't plotted by default.
    classes = sorted(igvDist, key=igvDist.get, reverse=True)
    total = sum(igvDist.values()) * 1.0
    writeCSV(os.path.join(outDir, sampleName + '_igv_dist_variant_level.csv'), "x,y\n", "{},{}\n",
             [(x, y) for x, y in zip(classes, map(lambda k: (igvDist[k] / total * 100), classes))])

    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igvDistSub = compressCountsGeneLevel(igvDist)
    #         for k in igvDist.keys():
    #             ksub = k.split('*')[0]
    #             igvDistSub[ksub] = igvDistSub.get(ksub, 0) + igvDist[k]
    plotDist(igvDistSub, sampleName, os.path.join(outDir, sampleName +
                                                  '_igv_dist_gene_level.csv'), rotateLabels=False, vertical=False,
             stream=stream)

    # Group IGVs based on the families and then write into a text file
    igvDistfam = compressCountsFamilyLevel(igvDistSub)
    #         for k in igvDistSub.keys():
    #             kfam = k.split('-')[0].split('/')[0]
    #             igvDistfam[kfam] = igvDistfam.get(kfam, 0) + igvDistSub[k]

    # Plot the family level distribution
    plotDist(igvDistfam, sampleName, os.path.join(outDir, sampleName +
                                                  '_igv_dist_family_level.csv'), stream=stream)

    # plot alignment length vs %identity
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'identity'],
                         ['Alignment Length', '%Identity'], os.path.join(outDir, sampleName +
                                                                         '_igv_align_quality_identity_hm.tsv'),
                         stream=stream)

    # plot alignment length vs bitScore
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'bitscore'],
                         ['Alignment Length', 'bitScore'], os.path.join(outDir, sampleName +
                                                                        '_igv_align_quality_bitscore_hm.tsv'),
                         stream=stream)

    # plot query start vs. subject start
    generateStatsHeatmap(stats, sampleName, ['vqstart', 'vstart'],
                         ['Query Start', 'Subject Start'], os.path.join(outDir, sampleName +
                                                                        '_igv_align_quality_start_hm.tsv'),
                         stream=stream)
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'vmismatches'],
                         ['Alignment Length', 'Mismatches'], os.path.join(outDir, sampleName +
                                                                          '_igv_align_quality_mismatches_hm.tsv'),
                         stream=stream)
    c = Counter(stats['vmismatches'].tolist())
    plotDist(c, sampleName, os.path.join(outDir, sampleName +
                                         '_igv_mismatches_dist.csv'), title='Number of Mismatches in V gene',
             proportion=True, rotateLabels=False, top=20, stream=stream)
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'vgaps'],
                         ['Alignment Length', 'Gaps'], os.path.join(outDir, sampleName +
                                                                    '_igv_align_quality_gaps_hm.tsv'), stream=stream)
    c = Counter(stats['vgaps'].tolist())
    plotDist(c, sampleName, os.path.join(outDir, sampleName +
                                         '_igv_gaps_dist.csv'), title='Number of Gaps in V gene',
             proportion=True, rotateLabels=False, top=20, stream=stream)
    #     print(np.percentile(stats, [0, 100], 0))
    #     summarizeStats(stats, outputDir+sampleName+'_stats_summary.txt')


def writeJAbundanceToFiles(stats, sampleName, outDir, stream=None):
    igjDist = Counter(stats["jgene"].tolist())
    igjDist = dict([(str(k), igjDist[k]) for k in igjDist])
    if len(igjDist) == 0:
        printto(stream, "WARNING: No IGJ hits were detected.", LEVEL.WARN)
        return

    plotDist(igjDist, sampleName, os.path.join(outDir, sampleName +
                                               '_igj_dist_variant_level.csv'), rotateLabels=False, vertical=False,
             stream=stream)

    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igjDistSub = compressCountsGeneLevel(igjDist)
    #     plotDist(igjDistSub, sampleName, outDir + sampleName +
    #              '_igj_dist_gene_level.csv', rotateLabels=False, vertical=False)
    #
    # Group IGVs based on the families and then write into a text file
    igjDistfam = compressCountsFamilyLevel(igjDistSub)
    # Plot the family level distribution
    plotDist(igjDistfam, sampleName, os.path.join(outDir, sampleName +
                                                  '_igj_dist_family_level.csv'),
             title='IGJ Abundance in Sample ' + sampleName, stream=stream)


def writeDAbundanceToFiles(stats, sampleName, outDir, stream=None):
    igdDist = Counter(stats["dgene"].tolist())
    igdDist = Counter(dict([(str(k), igdDist[k]) for k in igdDist]))
    if len(igdDist) == 0:
        printto(stream, "WARNING: No IGD hits were detected.", LEVEL.WARN)
        return

    # Write the counts of all IGVs into a text file
    # This isn't plotted by default, but we still write the csv file for it
    classes = sorted(igdDist, key=igdDist.get, reverse=True)
    total = sum(igdDist.values()) * 1.0
    writeCSV(os.path.join(outDir, sampleName + '_igd_dist_variant_level.csv'), "x,y\n", "{},{}\n",
             [(x, y) for x, y in zip(classes, map(lambda k: (igdDist[k] / total * 100), classes))])

    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igdDistSub = compressCountsGeneLevel(igdDist)
    plotDist(igdDistSub, sampleName, os.path.join(outDir, sampleName +
                                                  '_igd_dist_gene_level.csv'), rotateLabels=False, vertical=False,
             title='IGD Abundance in Sample ' + sampleName, stream=stream)

    # Group IGVs based on the families and then write into a text file
    igdDistfam = compressCountsFamilyLevel(igdDistSub)
    # Plot the family level distribution
    plotDist(igdDistfam, sampleName, os.path.join(outDir, sampleName +
                                                  '_igd_dist_family_level.csv'),
             title='IGD Abundance in Sample ' + sampleName, stream=stream)


def writeVJAssociationToFiles(stats, sampleName, outDir, stream=None):
    def canonicalFamilyName(v, j):
        vgene, jgene = str(v), str(j)
        if len(vgene) < 5 or len(jgene) < 5:
            return None, None
        return vgene.split("-")[0].split("/")[0].rstrip("D"), jgene.split("*")[0]

    fname = os.path.join(outDir, sampleName + "_vjassoc.csv")
    if os.path.exists(fname):
        return

    tally = defaultdict(lambda: defaultdict(int))
    for v, j in zip(stats['vgene'], stats['jgene']):
        vFamily, jFamily = canonicalFamilyName(v, j)
        if vFamily is None:  # jFamily is implicitly None too
            continue
        tally[vFamily][jFamily] += 1

    with open(fname, "w") as fp:
        writeBuffer = ""
        header = ["from", "to", "value"]
        writeBuffer += ",".join(header) + "\n"
        for vgene, dic in tally.items():
            for jgene, value in dic.items():
                writeBuffer += "{},{},{}\n".format(vgene, jgene, value)
        fp.write(writeBuffer)


def writeAbundanceToFiles(stats, sampleName, outDir, chain="hv", stream=None):
    writeVAbundanceToFiles(stats, sampleName, outDir, stream)
    writeJAbundanceToFiles(stats, sampleName, outDir, stream)
    writeVJAssociationToFiles(stats, sampleName, outDir, stream)
    if chain == "hv":
        writeDAbundanceToFiles(stats, sampleName, outDir, stream)
