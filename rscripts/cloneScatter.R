library(ggplot2)

options(scipen=999)
theme_set(theme_bw())

scatterPlot <- function(df1, df2, name1, name2, cloneClass) {
  df.union <- merge(df1, df2, by="Clonotype", all.y=TRUE, all.x=TRUE)
  # replace NaN with 0
  df.union[is.na(df.union)] <- 0
  
  # plot!
  gg <- ggplot(df.union, aes(x=Count.x, y=Count.y)) +
    geom_point() +
    labs(subtitle=paste(name2, " vs ", name1, " plot based on clonotype counts"),
         y=name2,
         x=name1,
         title=paste("Scatter plot of ", cloneClass, " clonotype counts"))
  return(gg)
}


scatterClones <- function(dataframes, sampleNames, outputPath, cloneClass) {
  # Plots a scatterplot of (i, i+1) samples. EG: given samples 1-3, plots
  # a scatterplot of pairs: (1,2) and (2,3).
  #
  # Args:
  #       dataframes: A list() type. List of dataframes
  #       sampleNames: A vector type. Vector of sample names (should correspond 1-1 with dataframes)
  #       outputPath: A string type. Path to save diagrams
  #       cloneClass: A string type. Title will specify what clonotype (e.g. CDR3/V-domain) this plot counts
  # Returns:
  #       Nothing. Plots n-1 scatter plots.
  nsamples <- length(dataframes)
  # this scatter plot doesn't make sense for 1 sample only
  stopifnot(nsamples > 1)
  stopifnot(nsamples == length(sampleNames))
  
  # dataframe cleanup
  colnames <- c("Clonotype", "Count")
  for (i in 1:(nsamples-1)) {
    df1 <- dataframes[[i]][colnames]
    df2 <- dataframes[[i+1]][colnames]
    p <- scatterPlot(df1, df2, sampleNames[i], sampleNames[i+1], cloneClass)
    ggsave(paste0(outputPath, sampleNames[i], "_vs_", sampleNames[i+1], "_clone_scatter.png"), plot = p, width = V_WIDTH, height = V_HEIGHT)
  }
}

#fs <- list.files(pattern = "PCR[123].*_cdr3_clonotypes_all_over.csv.gz$", recursive = TRUE, full.names = TRUE)
#scatterClones(lapply(fs, read.csv, stringsAsFactors = FALSE), c("PCR1_L001", "PCR2_L001", "PCR3_L001"), "/Users/harry/", "CDR3")