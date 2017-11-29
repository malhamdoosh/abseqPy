plotDist <- function(dataframes, sampleNames, plotTitle, vert = TRUE, xlabel = "", ylabel = "") {
  # Plots a dodging histogram for all sample in dataframes.
  #
  # Args:
  #   dataframes: A list() type. List of dataframes
  #   sampleNames: A vector type. Vector of sampleNames corresponding 1-1 to dataframes.
  #   plotTitle: A string type. 
  #   vert: A boolean type. True if the plot should be vertical
  #   xlabel: A string type. 
  #   ylabel: A string type.
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
    
    df <- dataframes[[i]]
    
    # for each dataframe, append a new column 'sample' to indicate where this df came from
    df$sample <- rep(sampleNames[i], nrow(df))
    
    # sort the values by descending order (strictly speaking,
    # not necessary because it's already sorted - but better be safe than sorry)
    if (vert) {
      dataframes[[i]] <- df[with(df, order(-y)), ]
    } else {
      dataframes[[i]] <- df[with(df, order(-x)), ]
    }
  }
  
  
  # make sure only the top CUTOFF is considered
  for (i in 1:frames) {
    df <- dataframes[[i]]
    if (nrow(df) > CUTOFF) {
      caps <- "Cutoff at top 15"
      dataframes[[i]] <- head(df, CUTOFF)
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
  
  if (missing(ylabel)) {
    ylabel <- "Proportion (%)"
  }
  
  if (!vert) {
    g <- ggplot(df.union, aes(y, x, label = sprintf("%0.2f%%", x))) + coord_flip()
    if (frames == 1) {
      g <- g + geom_text(hjust = -0.15, size = 3)
    }
  } else {
    g <- ggplot(df.union, aes(x, y, label = sprintf("%0.2f%%", y)))
    if (frames == 1) {
      g <- g + geom_text(vjust = -0.5, size = 3)
    }
    if (is.numeric(df.union$x)) {
      g <- g + scale_x_continuous(breaks = df.union$x, labels = df.union$x)
    }
  }
  g <- g + geom_bar(stat='identity', aes(fill = sample), width = 0.5, position = 'dodge') +
    #theme(text = element_text(size = 10), axis.text.x = element_text(angle = 74, vjust = 0.4)) +
    theme(text = element_text(size = 10)) +
    labs(title = plotTitle,
         x = xlabel,
         y = ylabel,
         caption = caps)
  return (g)
}
