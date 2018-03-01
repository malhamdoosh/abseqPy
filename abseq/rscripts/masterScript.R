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

# main functions are here
source(paste0(abSeqRoot, '/rscripts/annotAnalysis.R'))
source(paste0(abSeqRoot, '/rscripts/abundanceAnalysis.R'))
source(paste0(abSeqRoot, '/rscripts/productivityAnalysis.R'))
source(paste0(abSeqRoot, '/rscripts/diversityAnalysis.R'))
source(paste0(abSeqRoot, '/rscripts/primerSpecificity.R'))



# find all pairings (each different set of pairings is separated by a new line) -> reversed because we want to plot
# individual samples before doing sample comparisons
pairings <- rev(scan(metaFile, character(), quote = "", skip = 1))

# analysis types -> since ALL samples will have the SAME -t option, we can just infer from the first sample itself
analysis <- inferAnalyzed((unlist(strsplit(pairings[1], "\\?")[1]))[1])

# for each set of pairing - plot!
for (i in 1:length(pairings)) {
  ###                               Get information                                   ###
  
  # get pairings, split by "?", you'll have something like c("dirname1,dirname2,...", "sampleName1, sampleName2, ...")
  # i.e. pair[1] is a comma separated string of directories
  # i.e. pair[2] is a comma separated string of sample names
  pair <- unlist(strsplit(pairings[i], "\\?"))
  directories <- unlist(strsplit(pair[1], ","))
  sampleNames <- unlist(strsplit(pair[2], ","))
  
  # to get the result folder, we just need to look at any one of them, and take the full path until the penultimate directory
  decomposed <- unlist(strsplit(directories[1], "/"))
  resultFolder <- paste0(head(decomposed, n = length(decomposed) - 1), collapse = "/")
  
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
  if ('annot' %in% analysis) {
    annotOut <- paste0(outputDir, "annot/")
    if (!file.exists(annotOut)) {
      dir.create(annotOut)
    }
    annotDirectories <- unlist(lapply(directories, paste0, "/annot/"))
    annotAnalysis(annotDirectories, annotOut, sampleNames, mashedNames)   
  }
  
  ##################################################
  #                                                #
  #               ABUNDANCE PLOTS                  #
  #                                                #
  ##################################################
  if ('abundance' %in% analysis) {
    abunOut <- paste0(outputDir, "abundance/")
    if (!file.exists(abunOut)) {
      dir.create(abunOut)
    }
    abundanceDirectories <- unlist(lapply(directories, paste0, "/abundance/"))
    abundanceAnalysis(abundanceDirectories, abunOut, sampleNames, combinedNames, mashedNames)
  }
  
  ##################################################
  #                                                #
  #              PRODUCTIVITY PLOTS                #
  #                                                #
  ##################################################
  if ('productivity' %in% analysis) {
    prodOut <- paste0(outputDir, "productivity/")
    if (!file.exists(prodOut)) {
      dir.create(prodOut)
    }
    productivityDirectories <- unlist(lapply(directories, paste0, "/productivity/"))
    productivityAnalysis(productivityDirectories, prodOut, sampleNames, combinedNames, mashedNames)
  }
  
  ##################################################
  #                                                #
  #               DIVERSITY PLOTS                  #
  #                                                #
  ##################################################
  if ('diversity' %in% analysis) {
    diversityOut <- paste0(outputDir, "diversity/")
    if (!file.exists(diversityOut)) {
      dir.create(diversityOut)
    }
    diversityDirectories <- unlist(lapply(directories, paste0, "/diversity/"))
    diversityAnalysis(diversityDirectories, diversityOut, sampleNames, mashedNames)
  }
  
  ##################################################
  #                                                #
  #               PRIMER.S  PLOTS                  #
  #                                                #
  ##################################################
  if ('primer_specificity' %in% analysis) {
    primerOut <- paste0(outputDir, "primer_specificity/")
    if (!file.exists(primerOut)) {
      dir.create(primerOut)
    }
    primerDirectories <- unlist(lapply(directories, paste0, "/primer_specificity/"))
    args <- commandArgs(trailingOnly=TRUE)
    if (length(args) == 2) {
      primer5File <- args[1]
      primer3File <- args[2]
      primerAnalysis(primer5File, primer3File, primerDirectories, primerOut, combinedNames, mashedNames)
    }
  }
  
}

# print warnings
print(warnings())
