abundanceAnalysis <- function(abundanceDirectories, abunOut, sampleNames, combinedNames, mashedNames) {
  # where to find the files
  searchFiles <- listFilesInOrder(path = abundanceDirectories,
                     pattern = ".*ig[vdj]_dist_[family|gene|variant].*\\.csv(\\.gz)?$",
                     expectedRet = c(8, 5))            # 3 each from V and D, then 2 from J (no gene) or 5 (exclude the 3 from D)
  if (length(searchFiles) > 0) {
      abundancePlot(
        searchFiles,
        # what are the sample names (in-order)
        sampleNames,
        # output directory
        abunOut
      )
  }
  
  # plot igv mismatches distribution
  abunIgvMismatchFiles <- listFilesInOrder(path = abundanceDirectories, pattern = ".*_igv_mismatches_dist\\.csv(\\.gz)?$")
  if (length(abunIgvMismatchFiles) > 0) {
      subtitle <- paste("Total is", paste(lapply(abunIgvMismatchFiles, function(x) { as.integer(getTotal(x)) }), collapse = ", "))
      abunIgvMismatches <- plotDist(
        lapply(abunIgvMismatchFiles, read.csv, skip = 1),
        sampleNames,
        paste("Number of mismatches in V gene in", combinedNames),
        checkVert(abunIgvMismatchFiles[[1]]),
        subs = subtitle
      )
      ggsave(paste0(abunOut, mashedNames, "_igv_mismatches_dist.png"), plot = abunIgvMismatches, width = V_WIDTH, height = V_HEIGHT)
  }
  # plot igv gaps distribution
  abunIgvGapsFiles <- listFilesInOrder(path = abundanceDirectories, pattern = ".*_igv_gaps_dist\\.csv(\\.gz)?$")
  if (length(abunIgvGapsFiles) > 0) {
      subtitle <- paste("Total is", paste(lapply(abunIgvGapsFiles, function(x) { as.integer(getTotal(x)) }), collapse = ", "))
      abunIgvGaps <- plotDist(
        lapply(abunIgvGapsFiles, read.csv, skip = 1),
        sampleNames,
        paste("Number of gaps in V gene in ", combinedNames),
        checkVert(abunIgvGapsFiles[[1]]),
        subs = subtitle
      )
      ggsave(paste0(abunOut, mashedNames, "_igv_gaps_dist.png"), plot = abunIgvGaps, width = V_WIDTH, height = V_HEIGHT)
  }
  
  if (length(sampleNames) == 1) {
    # we can plot circlize if there's only one sample
    plotCirclize(sampleNames[1], abunOut)
  }
  
}