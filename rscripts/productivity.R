library(ggplot2)

productivityPlot <- function(dataframes, sampleNames) {
  # Plots a facet_grid histogram for all samples in dataframes.
  #
  # Args:
  #     dataframes: A list() type. List of dataframes
  #     sampleNames: A vector type. Vector of sampleNames corresponding 1-1 to dataframes
  # Returns:
  #     ggplot().
  
  nsamples <- length(dataframes)
  # sanity check
  stopifnot(length(sampleNames) == nsamples)
  
  #   ---- clean & pre-processing ----
  unusedCols <- c("X")
  # label each row with sample name and drop unused columns
  for (i in 1:nsamples) {
    df <- dataframes[[i]]
    df$round <- rep(sampleNames[i], nrow(df))
    dataframes[[i]] <- df[, !(names(df) %in% unusedCols)]
  }
  #   ---- done: clean & pre-processing ----
  
  # merge all samples
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
  g <- ggplot(df.union, aes(round, Percentage, label = sprintf("%0.2f%%", Percentage))) +
    geom_bar(stat="identity", aes(fill=Reason), width=0.5) +
    facet_grid(~ Productivity)+
    labs(title="Productivity",
         subtitle="Percentage of unproductive reads due to stop codons and frameshifts",
         x="Round",
         y="Percentage") +
    scale_y_continuous(limits=c(0,100)) +
    geom_text(position = position_stack(vjust = 0.5))
  return (g)
}
#files <- list.files(pattern = "PCR[123].*_productivity.csv$", full.names = TRUE, recursive = TRUE)
#dataframes <- lapply(files, read.csv, stringsAsFactors = FALSE)
#p <- productivityPlot(dataframes, c("PCR1_L001", "PCR2_L001", "PCR3_L001"))
#plot(p)
