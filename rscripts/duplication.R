library(ggplot2)
theme_set(theme_bw())

plotDuplication <- function(files, sampleNames, regions = c("CDR3", "V")) {
  # Plots duplication level
  # Args:
  #     files: A list() type. List of strings to _cdr_v_duplication.csv (pathname)
  #     sampleNames: A vector type. Vector of strings each representing the sample name
  #     regions: A vector type. Which regions to include in the plot. Default = c("CDR3", "V")
  # Returns:
  #     ggplot().
  nsamples <- length(files)
  # sanity check
  stopifnot(nsamples == length(sampleNames))
  
  # read xticks and xlimits from first 2 row
  ## little helper function - 'meta-function-composition' ##
  trimwsNoQuotes <- function(x) {
    gsub("'", "", trimws(x))
  }
  ## trimwsNoQuotes strips whitespaces AND single quotes ##
  fp <- file(files[[1]], "r")
  xticks <- strtoi(unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]], trimws)))
  xlabels <- unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]], trimwsNoQuotes))
  close(fp)
  
  # read files into dataframes
  dataframes <- lapply(files, read.csv, skip = 2)
  
  # pre-processing & cleanup
  for (i in 1:nsamples) {
    df <- dataframes[[i]]
    df$sample  <- rep(sampleNames[[i]], nrow(df))
    dataframes[[i]] <- df[df$region %in% regions, ]
  }
  
  # combine!
  if (nsamples > 1) {
    df.union <- rbind(dataframes[[1]], dataframes[[2]])
    if (nsamples > 2) {
      for (i in 3:nsamples) {
        df.union <- rbind(df.union, dataframes[[i]])
      }
    }
  } else {
    df.union <- dataframes[[1]]
  }
  
  g <- ggplot(df.union, aes(x = x, y = y)) + 
    geom_line(aes(linetype = region, color = sample)) + 
    scale_x_continuous(breaks = xticks, labels = xlabels) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Duplication levels of CDRs and V domains",
         x = "Duplication level",
         y = "Proportion of duplicated sequences")
}

# fs <- list.files(pattern = "PCR[123].*_cdr_v_duplication.csv", recursive = TRUE, full.names = TRUE)
# plot(plotDuplication(fs, c("PCR1", "PCR2", "PCR3"), c("CDR3", "CDR1", "CDR2", "V")))