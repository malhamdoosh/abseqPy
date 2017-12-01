# -------------------------------------------------------------------------------- #
#
#       Driver program - Plots all csv files with instructions taken from
#                        'metadata' file - rscripts_meta.tmp
#
# -------------------------------------------------------------------------------- #

# find AbSeq's root dir
metaFile <- "rscripts_meta.tmp"
f <- file(metaFile, "r")
abSeqRoot <- readLines(f, n = 1)
close(f)

# make sure all dependencies are met
source(paste0(abSeqRoot, '/rscripts/installer.R'))

library(ggplot2)

# source all required functions
source(paste0(abSeqRoot, '/rscripts/util.R'))
source(paste0(abSeqRoot, '/rscripts/plotDist.R'))
source(paste0(abSeqRoot, '/rscripts/summarySE.R'))
source(paste0(abSeqRoot, '/rscripts/abundance.R'))
source(paste0(abSeqRoot, '/rscripts/circlize.R'))
source(paste0(abSeqRoot, '/rscripts/cloneScatter.R'))
source(paste0(abSeqRoot, '/rscripts/clonotypeVenn.R'))
source(paste0(abSeqRoot, '/rscripts/duplication.R'))
source(paste0(abSeqRoot, '/rscripts/productivity.R'))
source(paste0(abSeqRoot, '/rscripts/rarefaction.R'))
source(paste0(abSeqRoot, '/rscripts/recapture.R'))
source(paste0(abSeqRoot, '/rscripts/regionAnalysis.R'))
source(paste0(abSeqRoot, '/rscripts/spectratype.R'))
source(paste0(abSeqRoot, '/rscripts/topNDist.R'))



# find all pairings (each different set of pairings is separated by a new line)
pairings <- scan(metaFile, character(), quote = "", skip = 1)

# for each set of pairing - plot!
for (i in 1:length(pairings)) {
  pair <- unlist(strsplit(pairings[i], "\\?"))
  directories <- unlist(strsplit(pair[1], ","))
  sampleNames <- unlist(strsplit(pair[2], ","))
  # resultFolder : this will be the same string as supplied / generated output folder in AbSeq's python code
  resultFolder <- unlist(strsplit(directories[1], "/"))[1]
  mashedNames <- paste(sampleNames, collapse = "_")
  combinedNames <- paste(sampleNames, collapse = ", ") 

  # different logic in obtaining folder names and sample directory names depending on sample lengths
  if (length(sampleNames) > 1) {
    outputDir <- paste0(resultFolder, '/', paste(sampleNames, collapse = "_vs_"), '/')
    dir.create(outputDir)
  } else {
    outputDir <- paste0(directories, '/')
  }
  
  # a little explaination: - the only variables you need to care about:
  
  # directories - a vector of directories for this pairing (could be a singleton vector) (root directory for each sample within a pairing)
  # eg: c("PCR1_BHZ123-AGGGA-ACGTA_L001", "PCR2_BHZ123-ACGT_L001", "PCR3_BHZ123_AGCT-GGACT_L001")
  
  # sampleNames - canonical names used by AbSeq (a vector of them, sample length as directories)
  # eg (for the above directories vector): c("PCR1_L001", "PCR2_L001", "PCR3_L001")
  
  # mashedNames - collapsing sampleNames with underscores for output png/pdfs
  
  # Plotting begins now - a little more walkthrough.
  # for each analysis, we:
  #   1) create <analysis>Out as the output directory for the specific analysis
  #   2) specialized directories to <analysis>Directories by pasting "/<analysis>/" to each root sample directory
  #   3) using the specialized directores, we use listFilesInOrder (list.files returns alphabetically, but our samples may be jumbled up when the
  #      provided ordering isn't alphabetical order (i.e. sampleNames isn't 1-1 with dataframes)) with regex to easily (and safely) find the right samples for the right
  #      function (within the analysis)
  #   4) plot!
  
  ##################################################
  #                                                #
  #               ANNOT PLOTS                      #
  #                                                #
  ##################################################
  annotOut <- paste0(outputDir, "annot/")
  dir.create(annotOut)
  annotDirectories <- unlist(lapply(directories, paste0, "/annot/"))
  
  # with outliers
  g <- plotSpectratype(
    lapply(listFilesInOrder(path = annotDirectories, pattern = ".*_all_clones_len_dist\\.csv(\\.gz)?$"), read.csv),
    sampleNames,
    title = "Sequence lengths",  # plotSpectratype names the sample(s) for us in the title; don't have to do it here
    xlabel = "Sequence Length (bp)",
    ylabel = "Proportion"
  )
  ggsave(paste0(annotOut, mashedNames, "_all_clones_len_dist.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  
  # without outliers
  g <- plotSpectratype(
    lapply(listFilesInOrder(path = annotDirectories, pattern = ".*_all_clones_len_dist_no_outliers\\.csv(\\.gz)?$"), read.csv),
    sampleNames,
    title = "Sequence lengths",  # plotSpectratype names the sample(s) for us in the title; don't have to do it here
    xlabel = "Sequence Length (bp)",
    ylabel = "Proportion"
  )
  ggsave(paste0(annotOut, mashedNames, "_all_clones_len_dist_no_outliers.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  
  ##################################################
  #                                                #
  #               ABUNDANCE PLOTS                  #
  #                                                #
  ##################################################
  abunOut <- paste0(outputDir, "abundance/")
  dir.create(abunOut)
  abundanceDirectories <- unlist(lapply(directories, paste0, "/abundance/"))
  abundancePlot(
    # where to find the files
    listFilesInOrder(path = abundanceDirectories,
                     pattern = ".*ig[vdj]_dist_[family|gene|variant].*\\.csv(\\.gz)?$",
                     expectedRet = 8),            # 3 each from V and D, then 2 from J (no gene)
    # what are the sample names (in-order)
    sampleNames,
    # output directory
    abunOut
  )
  
  # plot igv mismatches distribution
  abunIgvMismatchFiles <- listFilesInOrder(path = abundanceDirectories, pattern = ".*_igv_mismatches_dist\\.csv(\\.gz)?$")
  abunIgvMismatches <- plotDist(
    lapply(abunIgvMismatchFiles, read.csv, skip = 1),
    sampleNames,
    paste("Number of mismatches in V gene in", combinedNames),
    checkVert(abunIgvMismatchFiles[[1]])
  )
  ggsave(paste0(abunOut, mashedNames, "_igv_mismatches_dist.png"), plot = abunIgvMismatches, width = V_WIDTH, height = V_HEIGHT)
  # plot igv gaps distribution
  abunIgvGapsFiles <- listFilesInOrder(path = abundanceDirectories, pattern = ".*_igv_gaps_dist\\.csv(\\.gz)?$")
  abunIgvGaps <- plotDist(
    lapply(abunIgvGapsFiles, read.csv, skip = 1),
    sampleNames,
    paste("Number of gaps in V gene in ", combinedNames),
    checkVert(abunIgvGapsFiles[[1]])
  )
  ggsave(paste0(abunOut, mashedNames, "_igv_gaps_dist.png"), plot = abunIgvGaps, width = V_WIDTH, height = V_HEIGHT)
  
  if (length(sampleNames) == 1) {
    # we can plot circlize if there's only one sample
    plotCirclize(sampleNames[1], abunOut)
  }
  
  ##################################################
  #                                                #
  #              PRODUCTIVITY PLOTS                #
  #                                                #
  ##################################################
  prodOut <- paste0(outputDir, "productivity/")
  dir.create(prodOut)
  productivityDirectories <- unlist(lapply(directories, paste0, "/productivity/"))
  # main productivity file
  ########################
  prodFiles <- listFilesInOrder(path = productivityDirectories, pattern = ".*_productivity\\.csv(\\.gz)?$")
  g <- productivityPlot(lapply(prodFiles, read.csv, stringsAsFactors = FALSE), sampleNames)
  ggsave(paste0(prodOut, mashedNames, "_productivity.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  
  # sub-productivity files
  ########################
  regions <- c("cdr1", "cdr2", "cdr3", "fr1", "fr2", "fr3", "igv", "igd", "igj")
  
  # gaps_dist plots only
  gapPlots <- prodDistPlot(productivityDirectories, sampleNames, "Gaps in" ,"_gaps_dist\\.csv(\\.gz)?$", regions)
  i <- 1
  for (region in regions) {
    ggsave(paste0(prodOut, mashedNames, "_", region, "_gaps_dist.png"), plot = gapPlots[[i]], width = V_WIDTH, height = V_HEIGHT)
    i <- i + 1
  }
  
  # gaps_out_of_frame plots only (no igv, igd, ihj plots for this)
  subregions <- head(regions, n = 6)
  gapOutFramePlots <- prodDistPlot(productivityDirectories,
                                   sampleNames,
                                   "Gaps in",
                                   "_gaps_dist_out_of_frame\\.csv(\\.gz)?$",
                                   subregions)
  i <- 1
  for (region in subregions) {
    ggsave(paste0(prodOut, mashedNames, "_", region, "_gaps_dist_out_of_frame.png"), plot = gapOutFramePlots[[i]], width = V_WIDTH, height = V_HEIGHT)
    i <- i + 1
  }
  
  # mismatch dist only
  mismatchPlots <- prodDistPlot(productivityDirectories, sampleNames, "Mismatches in", "_mismatches_dist\\.csv(\\.gz)?$", regions)
  i <- 1
  for (region in regions) {
    ggsave(paste0(prodOut, mashedNames, "_", region, "_mismatches_dist.png"), plot = mismatchPlots[[i]], width = V_WIDTH, height = V_HEIGHT)
    i <- i + 1
  }
  
  # stop codon dist plot
  stopCodonPlot <- prodDistPlot(productivityDirectories, sampleNames,
                            "Stop codon in In-frame Clones",
                            "_stopcodon_dist_in_frame\\.csv(\\.gz)?$",
                            c("")) # no regions
  ggsave(paste0(prodOut, mashedNames, "_stopcodon_dist_in_frame.png"), plot = stopCodonPlot[[1]], width = V_WIDTH, height = V_HEIGHT)
  
  # vjframe plot
  vjframePlot <- prodDistPlot(productivityDirectories, sampleNames,
                              "V-D-J Rearrangement",
                              "_vjframe_dist\\.csv(\\.gz)?$",
                              c("")) # no regions
  ggsave(paste0(prodOut, mashedNames, "_vjframe_dist.png"), plot = vjframePlot[[1]], width = V_WIDTH, height = V_HEIGHT)
  
  # 3 special cases for IGV region
  # igv - inframe_unproductive, out of frame, productive
  igvInframeUnProductive <- prodDistPlot(productivityDirectories, sampleNames,
                                      "Abundance of In-frame Unproductive Clones in",
                                      "_dist_inframe_unproductive\\.csv(\\.gz)?$",
                                      c("igv"))
  ggsave(paste0(prodOut, mashedNames, "_igv_dist_inframe_unproductive.png"), plot = igvInframeUnProductive[[1]], width = V_WIDTH, height = V_HEIGHT)
  
  igvOutOfFrame <- prodDistPlot(productivityDirectories, sampleNames,
                                "Abundance of Out-Of-Frame Clones in",
                                "_dist_out_of_frame\\.csv(\\.gz)?$",
                                c("igv"))
  ggsave(paste0(prodOut, mashedNames, "_igv_dist_out_of_frame.png"), plot = igvOutOfFrame[[1]], width = V_WIDTH, height = V_HEIGHT)
  
  igvProductivity <- prodDistPlot(productivityDirectories, sampleNames,
                                  "Abundance of Productive Clones in",
                                  "_dist_productive\\.csv(\\.gz)?$",
                                  c("igv"))
  ggsave(paste0(prodOut, mashedNames, "_igv_dist_productive.png"), plot = igvProductivity[[1]], width = V_WIDTH, height = V_HEIGHT)
  
  ##################################################
  #                                                #
  #               CLONOTYPE PLOTS                  #
  #                                                #
  ##################################################
  diversityOut <- paste0(outputDir, "diversity/")
  dir.create(diversityOut)
  diversityDirectories <- unlist(lapply(directories, paste0, "/diversity/"))
  cdr3ClonesFile <- listFilesInOrder(path = diversityDirectories, pattern = ".*_cdr3_clonotypes_.*_over\\.csv(\\.gz)?$")
  if (length(sampleNames) > 1) {
    # plot scatter plot (CDR3 clonotypes)
    scatterClones(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
                  sampleNames,
                  diversityOut, "CDR3")
    
    # plot venn diagram (clonotypes)
    vennIntersection(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE), 
                     sampleNames,
                     paste0(diversityOut, mashedNames, "_cdr3_clonotypeIntersection.png"))
    
    # plot top N distribution (clonotypes)
    g <- topNDist(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
                  sampleNames)
    ggsave(paste0(diversityOut, mashedNames, "_top10Clonotypes.png"), plot = g, width = V_WIDTH_L, height = V_HEIGHT_L)
  }
  
  ##################################################
  #                                                #
  #               DIVERSITY PLOTS                  #
  #                                                #
  ##################################################
  # plot duplication, rarefaction, recapture
  for (reg in c("cdr", "cdr_v", "fr")) {
    if (reg == "cdr")
      includedRegions <- c("CDR1", "CDR2", "CDR3")
    else if (reg == "cdr_v")
      includedRegions <- c("CDR3", "V")
    else
      includedRegions <- c("FR1", "FR2", "FR3")
    
    # duplication
    g <- plotDuplication(listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_", reg,"_duplication\\.csv(\\.gz)?$")),
                         sampleNames,
                         includedRegions)
    ggsave(paste0(diversityOut, mashedNames, "_", reg, "_duplication.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
    
    # plot rarefaction
    g <- plotRarefaction(listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_", reg, "_rarefaction\\.csv(\\.gz)?$")),
                         sampleNames,
                         includedRegions)
    ggsave(paste0(diversityOut, mashedNames, "_", reg, "_rarefaction.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
    
    # plot recapture
    g <- plotRecapture(listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_", reg, "_recapture\\.csv(\\.gz)?$")),
                       sampleNames,
                       includedRegions)
    ggsave(paste0(diversityOut, mashedNames, "_", reg, "_recapture.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  }
  
  # MINISECTION:
  # ## SPECTRATYPES ###
  specOut <- paste0(diversityOut, "spectratypes/")
  dir.create(specOut)
  # CDR 1 - 3
  for (i in 1:3) {
    specFiles <- listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_cdr", i, "_spectratype\\.csv(\\.gz)?$"))
    g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                         sampleNames,
                         paste0("CDR", i))
    ggsave(paste0(specOut, mashedNames, "_cdr", i, "_spectratype.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  }
  
  # special case, no outliers plot for CDR3 only
  specFiles <- listFilesInOrder(path = diversityDirectories, pattern = ".*_cdr3_spectratype_no_outliers\\.csv(\\.gz)?$")
  g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                       sampleNames,
                       "CDR3")
  ggsave(paste0(specOut, mashedNames, "_cdr3_spectratype_no_outliers.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  
  # FR 1 - 4
  for (i in 1:4) {
    specFiles <- listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_fr", i, "_spectratype\\.csv(\\.gz)?$"))
    g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                         sampleNames,
                         paste0("FR", i))
    ggsave(paste0(specOut, mashedNames, "_fr", i, "_spectratype.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  }
  
  # entire V-domain
  specFiles <- listFilesInOrder(path = diversityDirectories, pattern = ".*_v_spectratype\\.csv(\\.gz)?$")
  g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                       sampleNames,
                       "V domain")
  ggsave(paste0(specOut, mashedNames, "_v_spectratype.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  
  # we can plot region analysis if there's only one sample
  #if (length(sampleNames) == 1) {
  #  # default = top 15
  #  g <- regionAnalysis(diversityOut, sampleNames[1])
  #  ggsave(paste0(diversityOut, mashedNames, "_region_analysis.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  #}
}

# print warnings
#warnings()