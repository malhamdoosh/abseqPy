library(ggplot2)
source('./summarySE.R')
theme_set(theme_bw())


plotRarefaction <- function(files, sampleNames, regions = c("CDR3", "V")) {
  # Plots rarefaction
  # Args:
  #      files: A list() type. A list of files consisting of path to samples
  #      sampleNames: A vector type. A vector of strings, each being the name of samples in files
  #      regions: A vector type. A vector of strings - regions to be included. Defaults to c("CDR3", "V")
  # Returns:
  #      ggplot().
  nsamples <- length(files)
  # sanity check
  stopifnot(length(sampleNames) == nsamples)
  
  # find the minimum (max)xtick value from all the samples
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
  
  # read files
  dataframes <- lapply(files, read.csv, skip = 1)
  
  # pre-processing & cleaning
  for (i in 1:nsamples) {
    df <- dataframes[[i]]
    df$sample <- rep(sampleNames[[i]], nrow(df))
    df <- df[df$region %in% regions, ]
    dataframes[[i]] <- summarySE(df, measurevar = 'y', groupvars = c('x', 'region', 'sample'))
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
  
  # plot!
  g <- ggplot(df.union, aes(x = x, y = y)) +
    geom_line(aes(linetype = region, color = sample)) +
    scale_x_continuous(breaks = xticks, limits = c(head(xticks, n = 1), tail(xticks, n = 1))) +
    geom_ribbon(aes(ymin = y - ci, ymax = y + ci, fill = region), alpha = 0.2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = "Rarefaction of CDRs and V Domains",
         subtitle = "Mean number of deduplicated sequences with 95% confidence interval",
         x = 'Sample size', y = "Number of deduplicated sequences")
  return (g)
}

#fs <- list.files(pattern = "PCR[123].*_cdr_v_rarefaction.csv.gz", recursive = TRUE, full.names = TRUE)
#plot(plotRarefaction(fs, c("PCR1_L001", "PCR2_L001", "PCR3_L001")))