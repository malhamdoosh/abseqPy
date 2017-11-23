library(ggplot2)

#TODO: (any summarySE.R too!)
sapply(Sys.glob('~/Documents/repo/abseq/rscripts/*.R'), source)

#TODO: remove rscipts_meta.tmp!
pairings <- scan("rscripts_meta.tmp", character(), quote = "")

# for each set of pairing - plot!
for (i in 1:length(pairings)) {
  pair <- unlist(strsplit(pairings[i], "\\?"))
  directories <- unlist(strsplit(pair[1], ","))
  sampleNames <- unlist(strsplit(pair[2], ","))
  
  # different logic in obtaining folder names and sample directory names depending on sample lengths
  if (length(sampleNames) > 1) {
    resultFolder <- head(tail(unlist(strsplit(directories[1], "/")), n = 2), n = 1)
    outputDir <- paste(resultFolder, '/', paste(sampleNames, collapse = "_vs_"), '/', sep = "")
    dir.create(outputDir)
  } else {
    resultFolder <- head(tail(unlist(strsplit(directories[1], "/")), n = 2), n = 1)
    sampleDirectory <- tail(tail(unlist(strsplit(directories[1], "/")), n = 2), n = 1) 
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
  # plot duplication
  g <- plotDuplication(list.files(path = directories,
                                  pattern = ".*_cdr_v_duplication\\.csv(\\.gz)?$",
                                  recursive = TRUE,
                                  full.names = TRUE),
                       sampleNames)
  ggsave(paste(diversityOut, "duplication.png", sep = ""), plot = g)
  
  # plot rarefaction
  g <- plotRarefaction(list.files(path = directories,
                                  pattern = ".*_cdr_v_rarefaction\\.csv(\\.gz)?$",
                                  recursive = TRUE,
                                  full.names = TRUE),
                       sampleNames)
  ggsave(paste(diversityOut, "rarefaction.png", sep = ""), plot = g)
  
  # plot recapture
  g <- plotRecapture(list.files(path = directories,
                                pattern = ".*_cdr_v_recapture\\.csv(\\.gz)?$",
                                recursive = TRUE,
                                full.names = TRUE),
                     sampleNames)
  ggsave(paste(diversityOut, "capture_recapture.png", sep = ""), plot = g)
  
  # plot (CDR3 only) spectratype
  specFiles <- list.files(path = directories,
                          pattern = ".*_cdr3_spectratype_no.*\\.csv(\\.gz)?$",
                          recursive = TRUE,
                          full.names = TRUE)
  g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                       sampleNames,
                       "CDR3")
  ggsave(paste(diversityOut, "cdr3_spectratype_no_outliers.png", sep = ""), plot = g)
}
