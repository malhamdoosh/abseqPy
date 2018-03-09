plotDist <- function(dataframes, sampleNames, plotTitle, vert = TRUE, xlabel = "", ylabel = "", perc = TRUE, subs = "", sortDist = TRUE) {
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
  #   sortDist: sorts plot by descending distribution order, True by default
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
    if (sortDist) {
      if (vert) {
        dataframes[[i]] <- df[with(df, order(-y)), ]
      } else {
        dataframes[[i]] <- df[with(df, order(-x)), ]
      }
    } else {
      dataframes[[i]] <- df
    }
    
    # canonicalize NaN, NA, and N/As to NA.
    dataframes[[i]][is.na(dataframes[[i]])] <- "NA"
    dataframes[[i]][dataframes[[i]] == "N/A"] <- "NA"
  }
  
  originals <- list()
  # make sure only the top CUTOFF is considered
  for (i in 1:frames) {
    df <- dataframes[[i]]
    originals[[i]] <- cbind(df)
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

  # -- complete merging -- #

  # need to compensate for topN in each dataframe
  # 1. gather unique x/y
  # 2. if x/y is in original table but not in top N table, append row to df.union
  # 3. otherwise, do nothing
    
  if (frames > 1) {
    if (vert) {
        topNx <- unique(df.union$x)
        for (x_ in topNx) {
            for (i in 1:frames) {
                if (x_ %in% originals[[i]]$x && !(x_ %in% dataframes[[i]]$x)) {
                    newRow <- originals[[i]][originals[[i]]$x == x_, ]
                    df.union <- rbind(df.union, newRow)
                }
            }
        }
    } else {
        topNy <- unique(df.union$y)
        for (y_ in topNy) {
            for (i in 1:frames) {
                if (y_ %in% originals[[i]]$y && !(y_ %in% dataframes[[i]]$y)) {
                    newRow <- originals[[i]][originals[[i]]$y == y_, ]
                    df.union <- rbind(df.union, newRow)
                }
            }
        }
    }
  }
  
  if (missing(ylabel)) {
    ylabel <- "Proportion"
    if (perc) {
      ylabel <- paste(ylabel, "(%)")
    }
  }
  if (!vert) {
    if (sortDist) {
      g <- ggplot(df.union, aes(x = reorder(y, x), y = x, label = sprintf(placeholder, x))) + coord_flip()
    } else {
      df.union$y <- factor(df.union$y, levels = unique(df.union$y))
      g <- ggplot(df.union, aes(x = y, y = x, label = sprintf(placeholder, x))) + coord_flip()
    }
    
    if (frames == 1) {
      # single sample -> blue colour plot
      g <- g + geom_text(hjust = 0.50, vjust = -0.5, size = 3, angle = -90) + 
        geom_bar(stat='identity', aes(fill = sample), position = 'dodge', fill = BLUEHEX, show.legend = FALSE)
    } else {
      # multiple samples -> multi-coloured plot
      g <- g + geom_bar(stat='identity', aes(fill = sample), position = 'dodge')
    }
    
  } else {
    if (sortDist) {
      g <- ggplot(df.union, aes(x = reorder(x, -y), y = y, label = sprintf(placeholder, y)))
    } else {
      df.union$x <- factor(df.union$x, levels = unique(df.union$x))
      g <- ggplot(df.union, aes(x = x, y = y, label = sprintf(placeholder, y)))
    }
    
    if (frames == 1) {
      # single sample -> blue colour plot
      g <- g + geom_text(vjust = -0.5, size = 3) + 
        geom_bar(stat='identity', aes(fill = sample), position = 'dodge', fill = BLUEHEX, show.legend = FALSE)
    } else {
      # multiple samples -> multi-coloured plot
      g <- g + geom_bar(stat='identity', aes(fill = sample), position = 'dodge')
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
