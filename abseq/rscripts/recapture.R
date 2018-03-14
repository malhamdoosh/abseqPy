library(ggplot2)
theme_set(theme_bw())



plotRecapture <- function(files, sampleNames, regions = c("CDR3", "V")) {
  # Plots recapture
  # Args:
  #     files: A list() type. List of _cdr_v_recapture.csv.gz files.
  #     sampleNames: A vector type. A vector of strings each representing the name of samples in files.
  #     regions: A vector type. A vector of strings - regions to be included in the plot. defaults to c("CDR3", "V")
  # Returns:
  #     ggplot()
  nsamples <- length(files)
  # sanity check
  stopifnot(nsamples == length(sampleNames))
  
  # find minimum (max)xticks value from each sample
  fp <- file(files[[1]], "r")
  xticks <- strtoi(unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]], trimws)))
  close(fp)
  
  # if there are more
  if (nsamples > 1) {
    for (i in 2:nsamples) {
      fp <- file(files[[i]], "r")
      candidate <- strtoi(unlist(lapply(strsplit(readLines(fp, n = 1), ",")[[1]], trimws)))
      close(fp)
      if (tail(candidate, n = 1) < tail(xticks, n = 1)) {
        xticks <- candidate
      }
    }
  }
  
  # read dataframes
  dataframes <- lapply(files, read.csv, skip = 1)
  
  # cleanup & pre-processing
  for (i in 1:nsamples) {
    df <- dataframes[[i]]
    # append sample name to a new column named sample
    df$sample <- rep(sampleNames[[i]], nrow(df))
    # only want selected regions - ignore others
    df <- df[df$region %in% regions,]
    # get mean, sd, se, and ci
    dataframes[[i]] <- summarySE(df, measurevar = 'y', groupvars = c("x", "region", "sample"))
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
  
  # make compound column for geom_ribbon (region . sample)
  df.union$compound <- paste(df.union$region, df.union$sample)
  
  # plot!
  p <- ggplot(df.union, aes(x = x, y = y))
  
  if (nsamples == 1) {
     p <- p + geom_line(aes(linetype = region, color = sample), color = BLUEHEX) + guides(color = FALSE)
  } else {
     p <- p + geom_line(aes(linetype = region, color = sample))
  }
  
  p <- p + scale_x_continuous(breaks = xticks) + 
    geom_ribbon(aes(ymin = y - ci, ymax = y + ci, fill = compound), alpha = 0.1, show.legend = FALSE) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste("Percent recapture of", paste(regions, collapse = ", "), "in", paste(sampleNames, collapse = ", ")),
         subtitle = "Mean number of recaptured sequences with 95% confidence interval",
         x = "Sample size", y = "Percent Recapture")
  return (p)
}

#fs <- list.files(pattern = "PCR[123].*_cdr_v_recapture.csv.gz$", recursive = TRUE, full.names = TRUE)
#plot(plotRecapture(fs, c("PCR1", "PCR2", "PCR3"), c("CDR1", "CDR2", "CDR3", "V")))