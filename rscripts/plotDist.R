plotDist <- function(dataframes, sampleNames, plotTitle, vert = TRUE, xlabel = "", ylabel = "", perc = TRUE, subs = "") {
  # Plots a dodging histogram for all sample in dataframes. If length(sampleNames) == 1, then the bars will
  # also have y-values (or x if horizontal plot) labels on them. Use 'perc' to control if the values are percentages.
  # Args:
  #   dataframes: A list() type. List of dataframes
  #   sampleNames: A vector type. Vector of sampleNames corresponding 1-1 to dataframes.
  #   plotTitle: A string type. 
  #   vert: A boolean type. True if the plot should be vertical
  #   xlabel: A string type. 
  #   ylabel: A string type.
  #   perc: Boolean type. True if data's axis is a % proportion (instead of 0-1)
  #         only meaningful if length(sampleNames) == 1
  #   subs: subtitle. defaults to empty
  # Returns:
  #   ggplot().
  
  frames <- length(dataframes)
  # sanity check
  stopifnot(frames == length(sampleNames))
  
  # define maximum cutoff x-axis values to display
  CUTOFF <- 15
  # If there was a cutoff, caps will display the right message, by default it's empty
  caps <- ""
  # if it was percentage, put a % sign in labels, otherwise just normal values
  if (perc) {
    placeholder <- "%0.2f%%"
  } else {
    placeholder <- "%0.2f"
  }


  
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
    ylabel <- "Proportion"
    if (perc) {
      ylabel <- paste(ylabel, "(%)")
    }
  }
  
  if (!vert) {
    g <- ggplot(df.union, aes(x = reorder(y, x), y = x, label = sprintf(placeholder, x))) + coord_flip()
    
    if (frames == 1) {
      # single sample -> blue colour plot
      g <- g + geom_text(hjust = -0.15, size = 3) + 
        geom_bar(stat='identity', aes(fill = sample), width = 0.5, position = 'dodge', fill = BLUEHEX, show.legend = FALSE)
    } else {
      # multiple samples -> multi-coloured plot
      g <- g + geom_bar(stat='identity', aes(fill = sample), width = 0.5, position = 'dodge')
    }
    
  } else {
    g <- ggplot(df.union, aes(x = reorder(x, -y), y = y, label = sprintf(placeholder, y)))
    
    if (frames == 1) {
      # single sample -> blue colour plot
      g <- g + geom_text(vjust = -0.5, size = 3) + 
        geom_bar(stat='identity', aes(fill = sample), width = 0.5, position = 'dodge', fill = BLUEHEX, show.legend = FALSE)
    } else {
      # multiple samples -> multi-coloured plot
      g <- g + geom_bar(stat='identity', aes(fill = sample), width = 0.5, position = 'dodge')
    }
  }
  g <- g + theme(text = element_text(size = 10)) +
    labs(title = plotTitle,
         subtitle = subs,
         x = xlabel,
         y = ylabel,
         caption = caps)
  return (g)
}
