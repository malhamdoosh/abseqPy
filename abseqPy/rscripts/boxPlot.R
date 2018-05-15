boxPlot <- function(dataframes, sampleNames, plotTitle, xlabel = "", ylabel = "", subs = "") {
  frames <- length(dataframes)
  
  stopifnot(frames == length(sampleNames))
  
  # add samplename into new "sample" column
  for (i in 1:frames) {
    dataframes[[i]]$sample <- rep(sampleNames[i], nrow(dataframes[[i]]))
  }
  
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
  if (frames == 1) {
    g <- ggplot(df.union, aes(x=x, y=y)) + geom_boxplot(varwidth = T, fill = BLUEHEX)
  } else {
    g <- ggplot(df.union, aes(x=sample, y=y)) + geom_boxplot(varwidth = T, fill = BLUEHEX) + facet_grid(~x) + theme(axis.text.x = element_text(angle = 75, hjust = 1))
  }
  
  g <- g + labs(title = plotTitle, sutitle = subs, x = xlabel, y = ylabel)
  return (g)
}