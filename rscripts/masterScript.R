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
  resultFolder <- unlist(strsplit(directories[1], "/"))[1]

  # different logic in obtaining folder names and sample directory names depending on sample lengths
  if (length(sampleNames) > 1) {
    outputDir <- paste(resultFolder, '/', paste(sampleNames, collapse = "_vs_"), '/', sep = "")
    dir.create(outputDir)
  } else {
    sampleDirectory <- unlist(strsplit(directories[1], "/"))[2]
    outputDir <- paste(resultFolder, '/', sampleDirectory, '/', sep = "")
  }
  
  ##################################################
  #                                                #
  #               ABUNDANCE PLOTS                  #
  #                                                #
  ##################################################
  abunOut <- paste(outputDir, "abundance/", sep = "")
  dir.create(abunOut)
  abundancePlot(
    # where to find the files
    list.files(path = directories,
               pattern = ".*ig[vdj]_dist_[family|gene|variant].*\\.csv(\\.gz)?$",
               full.names = TRUE,
               recursive = TRUE),
    # what are the sample names (in-order)
    sampleNames,
    # output directory
    abunOut
  )
  
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
  prodFiles <- list.files(path = directories,
                          pattern = ".*_productivity\\.csv(\\.gz)?$",
                          full.names = TRUE,
                          recursive = TRUE)
  g <- productivityPlot(lapply(prodFiles, read.csv, stringsAsFactors = FALSE),
                        sampleNames)
  ggsave(paste(prodOut, "productivity.png", sep = ""), plot = g)
  
  ##################################################
  #                                                #
  #               CLONOTYPE PLOTS                  #
  #                                                #
  ##################################################
  diversityOut <- paste(outputDir, "diversity/", sep = "")
  dir.create(diversityOut)
  cdr3ClonesFile <- list.files(path = directories,
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
                     paste(diversityOut, "cdr3_clonotypeIntersection.png", sep = ""))
    
    # plot top N distribution (clonotypes)
    g <- topNDist(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
                  sampleNames)
    ggsave(paste(diversityOut, "top10Clonotypes.png", sep = ""), plot = g)
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
    g <- plotDuplication(list.files(path = directories,
                                    pattern = paste(".*_", reg,"_duplication\\.csv(\\.gz)?$", sep = ""),
                                    recursive = TRUE,
                                    full.names = TRUE),
                         sampleNames, includedRegions)
    ggsave(paste(diversityOut, reg, "_duplication.png", sep = ""), plot = g)
    
    # plot rarefaction
    g <- plotRarefaction(list.files(path = directories,
                                    pattern = paste(".*_", reg, "_rarefaction\\.csv(\\.gz)?$", sep = ""),
                                    recursive = TRUE,
                                    full.names = TRUE),
                         sampleNames, includedRegions)
    ggsave(paste(diversityOut, reg, "_rarefaction.png", sep = ""), plot = g)
    
    # plot recapture
    g <- plotRecapture(list.files(path = directories,
                                  pattern = paste(".*_", reg, "_recapture\\.csv(\\.gz)?$", sep = ''),
                                  recursive = TRUE,
                                  full.names = TRUE),
                       sampleNames, includedRegions)
    ggsave(paste(diversityOut, reg, "_recapture.png", sep = ""), plot = g)
  }
  
  # MINISECTION:
  # ## SPECTRATYPES ###
  specOut <- paste(diversityOut, "spectratypes/", sep = "")
  dir.create(specOut)
  # CDR 1 - 3
  for (i in 1:3) {
    specFiles <- list.files(path = directories,
                            pattern = paste(".*_cdr", i, "_spectratype\\.csv(\\.gz)?$", sep = ""),
                            recursive = TRUE,
                            full.names = TRUE)
    g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                         sampleNames,
                         paste("CDR", i, sep = ""))
    ggsave(paste(specOut, "cdr", i, "_spectratype.png", sep = ""), plot = g)
  }
  
  # special case, no outliers plot for CDR3 only
  specFiles <- list.files(path = directories,
                          pattern = ".*_cdr3_spectratype_no.*\\.csv(\\.gz)?$",
                          recursive = TRUE,
                          full.names = TRUE)
  g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                       sampleNames,
                       "CDR3")
  ggsave(paste(specOut, "cdr3_spectratype_no_outliers.png", sep = ""), plot = g)
  
  # FR 1 - 4
  for (i in 1:4) {
    specFiles <- list.files(path = directories,
                            pattern = paste(".*_fr", i, "_spectratype\\.csv(\\.gz)?$", sep = ""),
                            recursive = TRUE,
                            full.names = TRUE)
    g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                         sampleNames,
                         paste("FR", i, sep = ""))
    ggsave(paste(specOut, "fr", i, "_spectratype.png", sep = ""), plot = g)
  }
  
  # entire V-domain
  specFiles <- list.files(path = directories,
                          pattern = paste(".*_v_spectratype\\.csv(\\.gz)?$", sep = ""),
                          recursive = TRUE,
                          full.names = TRUE)
  g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                       sampleNames,
                       "V domain")
  ggsave(paste(specOut, "v_spectratype.png", sep = ""), plot = g)
  
  # we can plot region analysis if there's only one sample
  if (length(sampleNames) == 1) {
    # default = top 15
    g <- regionAnalysis(diversityOut, sampleNames[1])
    ggsave(paste(diversityOut, "region_analysis.png", sep = ""), plot = g)
  }
}
