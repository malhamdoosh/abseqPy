library(ggplot2)


plotDist <- function(dataframes, sampleNames, gene, title="family") {
  # Plots a dodging histogram for all sample in dataframes.
  #
  # Args:
  #   dataframes: A list() type. List of dataframes
  #   sampleNames: A vector type. Vector of sampleNames corresponding 1-1 to dataframes.
  #   gene: A string type. The expression level of dataframe (either Family, Gene, Variant)
  #   title: A string type. If it's anything other than family, there will be a capped imposed on
  #          the number of x-axis values displayed (15 default)
  # Returns:
  #   ggplot().
  
  frames <- length(dataframes)
  # sanity check
  stopifnot(frames == length(sampleNames))
  
  # define maximum cutoff x-axis values to display
  CUTOFF <- 15
  # If there was a cutoff, caps will display the right message, by default it's empty
  caps <- ""
  
  # ------------  cleaning & transforming data   ----------- #
  
  for (i in 1:frames) {
    
    # remove the row that says 'TOTAL'
    df <- dataframes[[i]]
    df <- df[!(df$Germline.group == 'TOTAL'), ]
    
    # for each dataframe, append a new column 'sample' to indicate where this df came from
    df$sample <- rep(sampleNames[i], nrow(df))
    
    # sort the values by descending order (strictly speaking,
    # not necessary because it's already sorted - but better be safe than sorry)
    dataframes[[i]] <- df[with(df, order(-Count)), ]
  }
  
  
  # check if cut-off needs to happen (if title == family, there is always < 15 families anyway)
  if (title != "family") {
    for (i in 1:frames) {
      df <- dataframes[[i]]
      if (nrow(df) > CUTOFF) {
        caps <- "Cutoff at top 15"
        dataframes[[i]] <- head(df, CUTOFF)
      }
    }
  }
  
  # --------- complete: cleaning & transforming data -------- #
  
  
  # merge dataframes into one
  if (frames > 1) {
    df.union <- rbind(dataframes[[1]], dataframes[[2]])
    if (frames > 2) {
      for (i in 3:frames) {
        df <- dataframes[[i]]
        df.union <- rbind(df.union, df)
      }
    }
  } else {
    df.union <- dataframes[[1]]
  }
  
  g <- ggplot(df.union, aes(Germline.group, Percentage....)) +
    geom_bar(stat='identity', aes(fill = sample), width = 0.5, position = 'dodge') +
    theme(text = element_text(size = 10), axis.text.x = element_text(angle = 74, vjust = 0.4)) +
    labs(title = paste("Histogram of", gene, title, "distribution"),
         subtitle = paste("Percentage of IGH", gene, " ", title, " distribution across all samples", sep = ""),
         x = paste(gene, "-Germline ", title, sep = ""),
         y = "Percentage",
         caption = caps)
  return (g)
}


# sampleNames can be obtained as such : names(data_fames) <- gsub("\\.csv$", "", my_files)
abundancePlot <- function(fs, sampleNames) {
  for (expression in c("family", "gene")) {
    for (gene in c('v', 'd', 'j')) {
      # correction, J has no "gene" but rather variant
      if (gene == 'j' && expression == 'gene') {
        expression <- 'variant'
      }
      reg <- paste(".*_ig", gene, "_dist_", expression, "_level.csv$", sep = "")
      selectedFiles <- fs[grepl(reg, fs)]
      dataframes <- lapply(selectedFiles, read.csv, stringsAsFactors=FALSE)
      p <- plotDist(dataframes, sampleNames, toupper(gene), expression)
    }
  }
}

# https://stackoverflow.com/questions/17499013/how-do-i-make-a-list-of-data-frames
# test driver:
# fs can be obtained as such fs <- list.files(pattern = "\\.csv$")
#fs <- list.files(pattern = "PCR[123].*ig[vdj]_dist_[family|gene|variant].*\\.csv$", full.names = TRUE, recursive = TRUE)
#abundancePlot(fs, c("PCR1", "PCR2", "PCR3"))