'''
    Short description: Quality Control Analysis of Immunoglobulin Repertoire NGS (Paired-End MiSeq)    
    Author: Monther Alhamdoosh    
    Python Version: 2.7
    Changes log: check git commits. 
''' 

from collections import Counter, defaultdict
from IgRepReporting.igRepPlots import plotDist, generateStatsHeatmap
from IgRepertoire.igRepUtils import writeCountsToFile, compressCountsGeneLevel,\
    compressCountsFamilyLevel


def writeVAbundanceToFiles(stats, sampleName, outDir):
    igvDist = Counter(stats["vgene"].tolist())
    if (len(igvDist) == 0):
        print("WARNING: No IGV hits were detected.")
        return        
    
    # Write the counts of all IGVs into a text file
    writeCountsToFile(igvDist, outDir + sampleName + 
                      '_igv_dist_variant_level.csv')
    
    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igvDistSub = compressCountsGeneLevel(igvDist)
#         for k in igvDist.keys():
#             ksub = k.split('*')[0]
#             igvDistSub[ksub] = igvDistSub.get(ksub, 0) + igvDist[k]
    writeCountsToFile(igvDistSub, outDir + sampleName + 
                      '_igv_dist_gene_level.csv')
    plotDist(igvDistSub, sampleName, outDir + sampleName + 
             '_igv_dist_gene_level.png', rotateLabels=False, vertical=False)
    
    # Group IGVs based on the families and then write into a text file
    igvDistfam = compressCountsFamilyLevel(igvDistSub)
#         for k in igvDistSub.keys():
#             kfam = k.split('-')[0].split('/')[0]
#             igvDistfam[kfam] = igvDistfam.get(kfam, 0) + igvDistSub[k]
    writeCountsToFile(igvDistfam, outDir + sampleName + 
                       '_igv_dist_family_level.csv')
    
    # Plot the family level distribution
    plotDist(igvDistfam, sampleName, outDir + sampleName + 
             '_igv_dist_family_level.png')
    
    # plot alignment length vs %identity
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'identity'],
                     ['Alignment Length', '%Identity'] , outDir + sampleName + 
              '_igv_align_quality_identity_hm.png')
    # plot alignment length vs bitScore
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'bitscore'],
                     ['Alignment Length', 'bitScore'] , outDir + sampleName + 
              '_igv_align_quality_bitscore_hm.png')
    # plot query start vs. subject start
    generateStatsHeatmap(stats, sampleName, ['vqstart', 'vstart'],
                     ['Query Start', 'Subject Start'] , outDir + sampleName + 
              '_igv_align_quality_start_hm.png')
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'vmismatches'],
                     ['Alignment Length', 'Mismatches'] , outDir + sampleName + 
              '_igv_align_quality_mismatches_hm.png')
    c = Counter(stats['vmismatches'].tolist())
    plotDist(c, sampleName, outDir + sampleName + 
             '_igv_mismatches_dist.png', title='Number of Mismatches in V gene',
             proportion=True, rotateLabels=False, top=20) 
    generateStatsHeatmap(stats, sampleName, ['alignlen', 'vgaps'],
                     ['Alignment Length', 'Gaps'] , outDir + sampleName + 
              '_igv_align_quality_gaps_hm.png')
    c = Counter(stats['vgaps'].tolist())
    plotDist(c, sampleName, outDir + sampleName + 
             '_igv_gaps_dist.png', title='Number of Gaps in V gene',
             proportion=True, rotateLabels=False, top=20) 
    #     print(np.percentile(stats, [0, 100], 0))
    #     summarizeStats(stats, outputDir+sampleName+'_stats_summary.txt')


def writeJAbundanceToFiles(stats, sampleName, outDir):
    igjDist = Counter(stats["jgene"].tolist())
    igjDist = {str(k) : igjDist[k] for k in igjDist}
    if (len(igjDist) == 0):
        print("WARNING: No IGJ hits were detected.")
        return        
    
    # Write the counts of all IGVs into a text file
    writeCountsToFile(igjDist, outDir + sampleName + 
                      '_igj_dist_variant_level.csv')
    plotDist(igjDist, sampleName, outDir + sampleName + 
                  '_igj_dist_variant_level.png', rotateLabels=False, vertical=False)
    
    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igjDistSub = compressCountsGeneLevel(igjDist)
#     writeCountsToFile(igjDistSub, outDir + sampleName + 
#                       '_igj_dist_gene_level.csv')
#     plotDist(igjDistSub, sampleName, outDir + sampleName + 
#              '_igj_dist_gene_level.png', rotateLabels=False, vertical=False)
#     
    # Group IGVs based on the families and then write into a text file
    igjDistfam = compressCountsFamilyLevel(igjDistSub)
    writeCountsToFile(igjDistfam, outDir + sampleName + 
                       '_igj_dist_family_level.csv')    
    # Plot the family level distribution
    plotDist(igjDistfam, sampleName, outDir + sampleName + 
             '_igj_dist_family_level.png',
             title = 'IGJ Abundance in Sample ' + sampleName )


def writeDAbundanceToFiles(stats, sampleName, outDir):
    igdDist = Counter(stats["dgene"].tolist())
    igdDist = Counter({str(k) : igdDist[k] for k in igdDist})
    if (len(igdDist) == 0):
        print("WARNING: No IGD hits were detected.")
        return        
    
    # Write the counts of all IGVs into a text file
    writeCountsToFile(igdDist, outDir + sampleName + 
                      '_igd_dist_variant_level.csv')
    
    # Group IGVs based on the subfamilies (gene level) and then write into a text file
    igdDistSub = compressCountsGeneLevel(igdDist)
    writeCountsToFile(igdDistSub, outDir + sampleName + 
                      '_igd_dist_gene_level.csv')
    plotDist(igdDistSub, sampleName, outDir + sampleName + 
             '_igd_dist_gene_level.png', rotateLabels=False, vertical=False,
             title = 'IGD Abundance in Sample ' + sampleName )
    
    # Group IGVs based on the families and then write into a text file
    igdDistfam = compressCountsFamilyLevel(igdDistSub)
    writeCountsToFile(igdDistfam, outDir + sampleName + 
                       '_igd_dist_family_level.csv')    
    # Plot the family level distribution
    plotDist(igdDistfam, sampleName, outDir + sampleName + 
             '_igd_dist_family_level.png',
             title = 'IGD Abundance in Sample ' + sampleName)


def writeVJAssociationToFiles(stats, sampleName, outDir):
    def canonicalFamilyName(v, j):
        vgene, jgene = str(v), str(j)
        if len(vgene) < 5 or len(jgene) < 5:
            return None, None
        return vgene.split("-")[0].split("/")[0].rstrip("D"), jgene.split("*")[0]

    tally = defaultdict(lambda: defaultdict(int))
    for v, j in zip(stats['vgene'], stats['jgene']):
        vFamily, jFamily = canonicalFamilyName(v, j)
        if vFamily is None:     # jFamily is implicitly None too
            continue
        tally[vFamily][jFamily] += 1

    with open(outDir + sampleName + "_vjassoc.csv", "w") as fp:
        writeBuffer = ""
        header = ["from", "to", "value"]
        writeBuffer += ",".join(header) + "\n"
        for vgene, dic in tally.items():
            for jgene, value in dic.items():
                writeBuffer += "{},{},{}\n".format(vgene, jgene, value)
        fp.write(writeBuffer)


def writeAbundanceToFiles(stats, sampleName, outDir, chain = "hv"):
    writeVAbundanceToFiles(stats, sampleName, outDir)
    writeJAbundanceToFiles(stats, sampleName, outDir)
    writeVJAssociationToFiles(stats, sampleName, outDir)
    if (chain == "hv"):
        writeDAbundanceToFiles(stats, sampleName, outDir)
