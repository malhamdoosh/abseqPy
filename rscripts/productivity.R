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
    labs(title="Productivity proportions",
         subtitle="Percentage of unproductive reads due to stop codons and frameshifts",
         x="Round",
         y="Percentage (%)") +
    scale_y_continuous(limits=c(0,100)) +
    geom_text(position = position_stack(vjust = 0.5))
  return (g)
}

prodDistPlot <- function(productivityDirectories, sampleNames, title, reg, saveNames,
                         regions = c("cdr1", "cdr2", "cdr3", "fr1", "fr2", "fr3", "igv", "igd", "igj")) {
  # Plots a distribution plot for different productivity analysis files
  # Args:
  #     productivityDirectories: A vector type. directories where all productivity csv files lives (usually <samplename>/productivity/)
  #     sampleNames: A vector type. Vector of strings.
  #     title: Title to give the plot
  #     reg: Regular expression to find the right files for this particular distribution plot (see samples in masterScript.R)
  #     regions: Most of the dist plots are regional based. use c("") if no regions are involved
  # Returns:
  #     a list of ggplot()s.
  stopifnot(length(regions) == length(saveNames))
  i = 1
  for (region in regions) {
    fs <- listFilesInOrder(path = productivityDirectories, pattern = paste0(".*", region, reg))
    if (length(fs) > 0) {
        dataframes <- lapply(fs, read.csv, stringsAsFactors = FALSE, skip = 1)
        plotTitle <- paste(title, toupper(region), "in", paste(sampleNames, collapse = ", "))
        subtitle <- paste("Total is" , paste(lapply(fs, function(x) { as.integer(getTotal(x)) } ), collapse = ", "))
        g <- plotDist(dataframes, sampleNames, plotTitle, checkVert(fs[[1]]), subs = subtitle)
        ggsave(saveNames[i], plot = g, width = V_WIDTH, height = V_HEIGHT)
    }
    i <- i + 1
  }
}

#files <- list.files(pattern = "PCR[123].*_productivity.csv$", full.names = TRUE, recursive = TRUE)
#dataframes <- lapply(files, read.csv, stringsAsFactors = FALSE)
#p <- productivityPlot(dataframes, c("PCR1_L001", "PCR2_L001", "PCR3_L001"))
#plot(p)
