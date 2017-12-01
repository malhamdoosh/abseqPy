topNDist <- function(dataframes, sampleNames, top = 10) {
  # Plots a histogram of top N clonotypes
  # Args:
  #      dataframes: A list() type. List of dataframes.
  #      sampleNames: A vector type. vector of strings representing sample names
  #                   (should have one-to-one correspondence with dataframes)
  #      top: Top N clonotypes to plot
  # Returns: ggplot().
  nsamples <- length(dataframes)
  # sanity check
  stopifnot(nsamples == length(sampleNames))
  
  # --- cleanup & pre-processing ---
  colNames <- c("Clonotype", "Count")
  for (i in 1:nsamples) {
    df <- dataframes[[i]]
    # only need top N, also take the columns we're interested in only
    df <- head(df[colNames], top)
    # append sample name to distinguish data when merged later on
    df$round <- rep(sampleNames[i], nrow(df))
    # normalize percentage to top N
    df$Count <- df$Count / sum(df$Count)
    dataframes[[i]] <- df
  }
  
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
  # --- done: cleanup & pre-processing ---
  
  # plot!
  # colour suggestion taken from
  # https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes/41230685 and modified
  c30 <- c("dodgerblue2","#E31A1C", # red
           "green4",
           "#6A3D9A", # purple
           "#FF7F00", # orange
           "black","gold1",
           "skyblue2","#FB9A99", # lt pink
           "olivedrab1",
           "#CAB2D6", # lt purple
           "#FDBF6F", # lt orange
           "gray70", "khaki2",
           "maroon","orchid1","deeppink1","blue1","steelblue4",
           "darkturquoise","green1","yellow4","yellow3",
           "aquamarine", "darkorange4", "mediumpurple1", "dimgrey", "darkseagreen1", "lightyellow", "coral2")
  
  g <- ggplot(df.union, aes(x=round, y=Count)) +
    geom_bar(stat='identity', aes(fill=Clonotype)) +
    theme(legend.position="bottom", legend.box = "horizontal", legend.title=element_blank(), legend.text=element_text(size=7)) +
    labs(title=paste("Top", top, "clonotype across each sample"),
         subtitle=paste0("Colour coded clonotypes, distribution of each clonotype is relative to top ", top, ", not overall."),
         x="round",
         y="Distribution")  +
    scale_fill_manual(values=c30)
  return (g)
}

#fs <- list.files(pattern = "PCR[123].*_cdr3_clonotypes_full_over.csv.gz", recursive = TRUE, full.names = TRUE)
#p <- topNDist(lapply(fs, read.csv, stringsAsFactors = FALSE), c("PCR1_L001", "PCR2_L001", "PCR3_L001"))
#plot(p)