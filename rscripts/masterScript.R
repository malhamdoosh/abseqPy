# -------------------------------------------------------------------------------- #
#
#       Driver program - Plots all csv files with instructions taken from
#                        'metadata' file - rscripts_meta.tmp
#
# -------------------------------------------------------------------------------- #

library(ggplot2)

metaFile <- "rscripts_meta.tmp"
f <- file(metaFile, "r")
abSeqRoot <- readLines(f, n = 1)
close(f)

source(paste(abSeqRoot, '/rscripts/util.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/plotDist.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/summarySE.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/abundance.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/circlize.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/cloneScatter.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/clonotypeVenn.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/duplication.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/productivity.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/rarefaction.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/recapture.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/regionAnalysis.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/spectratype.R', sep = ""))
source(paste(abSeqRoot, '/rscripts/topNDist.R', sep = ""))



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

  # different logic in obtaining folder names and sample directory names depending on sample lengths
  if (length(sampleNames) > 1) {
    outputDir <- paste(resultFolder, '/', paste(sampleNames, collapse = "_vs_"), '/', sep = "")
    dir.create(outputDir)
  } else {
    outputDir <- paste(directories, '/', sep = "")
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
  #   3) using the specialized directores, we use list.files with regex to easily (and safely) find the right samples for the right
  #      function (within the analysis)
  #   4) plot!
  
  ##################################################
  #                                                #
  #               ANNOT PLOTS                      #
  #                                                #
  ##################################################
  annotOut <- paste(outputDir, "annot/", sep = "")
  dir.create(annotOut)
  annotDirectories <- unlist(lapply(directories, paste, "/annot/", sep = ""))
  g <- plotSpectratype(
    lapply(list.files(path = annotDirectories,
                      pattern = ".*_all_clones_len_dist\\.csv(\\.gz)?$",
                      full.names = TRUE,
                      recursive = TRUE),
           read.csv),
    sampleNames,
    title = "Sequence lengths",
    xlabel = "Sequence Length (bp)",
    ylabel = "Proportion"
  )
  ggsave(paste(annotOut, mashedNames, "_all_clones_len_dist.png", sep = ""), plot = g)
  
  
  ##################################################
  #                                                #
  #               ABUNDANCE PLOTS                  #
  #                                                #
  ##################################################
  abunOut <- paste(outputDir, "abundance/", sep = "")
  dir.create(abunOut)
  abundanceDirectories <- unlist(lapply(directories, paste, "/abundance/", sep = ""))
  abundancePlot(
    # where to find the files
    list.files(path = abundanceDirectories,
               pattern = ".*ig[vdj]_dist_[family|gene|variant].*\\.csv(\\.gz)?$",
               full.names = TRUE,
               recursive = TRUE),
    # what are the sample names (in-order)
    sampleNames,
    # output directory
    abunOut
  )
  
  # plotDist <- function(dataframes, sampleNames, plotTitle, vert = TRUE, xlabel = "", ylabel = "")
  # plot igv mismatches distribution
  abunIgvMismatchFiles <- list.files(path = abundanceDirectories,
                                     pattern = ".*_igv_mismatches_dist\\.csv(\\.gz)?$",
                                     full.names = TRUE,
                                     recursive = TRUE)
  abunIgvMismatches <- plotDist(
    lapply(abunIgvMismatchFiles, read.csv, skip = 1),
    sampleNames,
    "Number of mismatches in V gene",
    checkVert(abunIgvMismatchFiles[[1]])
  )
  ggsave(paste(abunOut, mashedNames, "_igv_mismatches_dist.png", sep = ""), plot = abunIgvMismatches)
  # plot igv gaps distribution
  abunIgvGapsFiles <- list.files(path = abundanceDirectories,
                                 pattern = ".*_igv_gaps_dist\\.csv(\\.gz)?$",
                                 full.names = TRUE,
                                 recursive = TRUE)
  abunIgvGaps <- plotDist(
    lapply(abunIgvGapsFiles, read.csv, skip = 1),
    sampleNames,
    "Number of gaps in V gene",
    checkVert(abunIgvGapsFiles[[1]])
  )
  ggsave(paste(abunOut, mashedNames, "_igv_gaps_dist.png", sep = ""), plot = abunIgvGaps)
  
  if (length(sampleNames) == 1) {
    # we can plot circlize if there's only one sample
    plotCirclize(sampleNames[1], abunOut)
  }
  
  ##################################################
  #                                                #
  #              PRODUCTIVITY PLOTS                #
  #                                                #
  ##################################################
  prodOut <- paste(outputDir, "productivity/", sep = "")
  dir.create(prodOut)
  productivityDirectories <- unlist(lapply(directories, paste, "/productivity/", sep = ""))
  # main productivity file
  ########################
  prodFiles <- list.files(path = productivityDirectories,
                          pattern = ".*_productivity\\.csv(\\.gz)?$",
                          full.names = TRUE,
                          recursive = TRUE)
  g <- productivityPlot(lapply(prodFiles, read.csv, stringsAsFactors = FALSE),
                        sampleNames)
  ggsave(paste(prodOut, mashedNames, "_productivity.png", sep = ""), plot = g)
  
  # sub-productivity files
  ########################
  regions <- c("cdr1", "cdr2", "cdr3", "fr1", "fr2", "fr3", "igv", "igd", "igj")
  
  # gaps_dist plots only
  gapPlots <- prodDistPlot(productivityDirectories, sampleNames, "Gaps in" ,"_gaps_dist\\.csv(\\.gz)?$", regions)
  i <- 1
  for (region in regions) {
    ggsave(paste(prodOut, mashedNames, "_", region, "_gaps_dist.png", sep = ""), plot = gapPlots[[i]])  
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
    ggsave(paste(prodOut, mashedNames, "_", region, "_gaps_dist_out_of_frame.png", sep = ""), plot = gapOutFramePlots[[i]])  
    i <- i + 1
  }
  
  # mismatch dist only
  mismatchPlots <- prodDistPlot(productivityDirectories, sampleNames, "Mismatches in", "_mismatches_dist\\.csv(\\.gz)?$", regions)
  i <- 1
  for (region in regions) {
    ggsave(paste(prodOut, mashedNames, "_", region, "_mismatches_dist.png", sep = ""), plot = mismatchPlots[[i]])  
    i <- i + 1
  }
  
  # stop codon dist plot
  stopCodonPlot <- prodDistPlot(productivityDirectories, sampleNames,
                            "Stop codon in In-frame Clones",
                            "_stopcodon_dist_in_frame\\.csv(\\.gz)?$",
                            c("")) # no regions
  ggsave(paste(prodOut, mashedNames, "_stopcodon_dist_in_frame.png", sep = ""), plot = stopCodonPlot[[1]])
  
  # vjframe plot
  vjframePlot <- prodDistPlot(productivityDirectories, sampleNames,
                              "V-D-J Rearrangement",
                              "_vjframe_dist\\.csv(\\.gz)?$",
                              c("")) # no regions
  ggsave(paste(prodOut, mashedNames, "_vjframe_dist.png", sep = ""), plot = vjframePlot[[1]])
  
  # 3 special cases for IGV region
  # igv - inframe_unproductive, out of frame, productive
  igvInframeUnProductive <- prodDistPlot(productivityDirectories, sampleNames,
                                      "Abundance of In-frame Unproductive Clones in",
                                      "_dist_inframe_unproductive\\.csv(\\.gz)?$",
                                      c("igv"))
  ggsave(paste(prodOut, mashedNames, "_igv_dist_inframe_unproductive.png", sep = ""), plot = igvInframeUnProductive[[1]])
  
  igvOutOfFrame <- prodDistPlot(productivityDirectories, sampleNames,
                                "Abundance of Out-Of-Frame Clones in",
                                "_dist_out_of_frame\\.csv(\\.gz)?$",
                                c("igv"))
  ggsave(paste(prodOut, mashedNames, "_igv_dist_out_of_frame.png", sep = ""), plot = igvOutOfFrame[[1]])
  
  igvProductivity <- prodDistPlot(productivityDirectories, sampleNames,
                                  "Abundance of Productive Clones in",
                                  "_dist_productive\\.csv(\\.gz)?$",
                                  c("igv"))
  ggsave(paste(prodOut, mashedNames, "_igv_dist_productive.png", sep = ""), plot = igvProductivity[[1]])
  
  ##################################################
  #                                                #
  #               CLONOTYPE PLOTS                  #
  #                                                #
  ##################################################
  diversityOut <- paste(outputDir, "diversity/", sep = "")
  dir.create(diversityOut)
  diversityDirectories <- unlist(lapply(directories, paste, "/diversity/", sep = ""))
  cdr3ClonesFile <- list.files(path = diversityDirectories,
                               pattern = ".*_cdr3_clonotypes_.*_over\\.csv(\\.gz)?$",
                               recursive = TRUE,
                               full.names = TRUE)
  if (length(sampleNames) > 1) {
    # plot scatter plot (CDR3 clonotypes)
    scatterClones(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
                  sampleNames,
                  diversityOut, "CDR3")
    
    # plot venn diagram (clonotypes)
    vennIntersection(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE), 
                     sampleNames,
                     paste(diversityOut, mashedNames, "_cdr3_clonotypeIntersection.png", sep = ""))
    
    # plot top N distribution (clonotypes)
    g <- topNDist(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
                  sampleNames)
    ggsave(paste(diversityOut, mashedNames, "_top10Clonotypes.png", sep = ""), plot = g)
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
    g <- plotDuplication(list.files(path = diversityDirectories,
                                    pattern = paste(".*_", reg,"_duplication\\.csv(\\.gz)?$", sep = ""),
                                    recursive = TRUE,
                                    full.names = TRUE),
                         sampleNames, includedRegions)
    ggsave(paste(diversityOut, mashedNames, "_", reg, "_duplication.png", sep = ""), plot = g)
    
    # plot rarefaction
    g <- plotRarefaction(list.files(path = diversityDirectories,
                                    pattern = paste(".*_", reg, "_rarefaction\\.csv(\\.gz)?$", sep = ""),
                                    recursive = TRUE,
                                    full.names = TRUE),
                         sampleNames, includedRegions)
    ggsave(paste(diversityOut, mashedNames, "_", reg, "_rarefaction.png", sep = ""), plot = g)
    
    # plot recapture
    g <- plotRecapture(list.files(path = diversityDirectories,
                                  pattern = paste(".*_", reg, "_recapture\\.csv(\\.gz)?$", sep = ''),
                                  recursive = TRUE,
                                  full.names = TRUE),
                       sampleNames, includedRegions)
    ggsave(paste(diversityOut, mashedNames, "_", reg, "_recapture.png", sep = ""), plot = g)
  }
  
  # MINISECTION:
  # ## SPECTRATYPES ###
  specOut <- paste(diversityOut, "spectratypes/", sep = "")
  dir.create(specOut)
  # CDR 1 - 3
  for (i in 1:3) {
    specFiles <- list.files(path = diversityDirectories,
                            pattern = paste(".*_cdr", i, "_spectratype\\.csv(\\.gz)?$", sep = ""),
                            recursive = TRUE,
                            full.names = TRUE)
    g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                         sampleNames,
                         paste("CDR", i, sep = ""))
    ggsave(paste(specOut, mashedNames, "_cdr", i, "_spectratype.png", sep = ""), plot = g)
  }
  
  # special case, no outliers plot for CDR3 only
  specFiles <- list.files(path = diversityDirectories,
                          pattern = ".*_cdr3_spectratype_no_outliers\\.csv(\\.gz)?$",
                          recursive = TRUE,
                          full.names = TRUE)
  g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                       sampleNames,
                       "CDR3")
  ggsave(paste(specOut, mashedNames, "_cdr3_spectratype_no_outliers.png", sep = ""), plot = g)
  
  # FR 1 - 4
  for (i in 1:4) {
    specFiles <- list.files(path = diversityDirectories,
                            pattern = paste(".*_fr", i, "_spectratype\\.csv(\\.gz)?$", sep = ""),
                            recursive = TRUE,
                            full.names = TRUE)
    g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                         sampleNames,
                         paste("FR", i, sep = ""))
    ggsave(paste(specOut, mashedNames, "_fr", i, "_spectratype.png", sep = ""), plot = g)
  }
  
  # entire V-domain
  specFiles <- list.files(path = diversityDirectories,
                          pattern = paste(".*_v_spectratype\\.csv(\\.gz)?$", sep = ""),
                          recursive = TRUE,
                          full.names = TRUE)
  g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                       sampleNames,
                       "V domain")
  ggsave(paste(specOut, mashedNames, "_v_spectratype.png", sep = ""), plot = g)
  
  # we can plot region analysis if there's only one sample
  #if (length(sampleNames) == 1) {
  #  # default = top 15
  #  g <- regionAnalysis(diversityOut, sampleNames[1])
  #  ggsave(paste(diversityOut, mashedNames, "_region_analysis.png", sep = ""), plot = g)
  #}
}

# print warnings
#warnings()