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
  palette <- c("dodgerblue2","#E31A1C", # red
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
  if (length(unique(df.union$Clonotype)) > 30) {
    print("WARNING: Too many unique clonotypes are being plotted - extrapolating palette for top10 clonotype dist plot.")
    getPalatte <- colorRampPalette(brewer.pal(8, 'Accent'))
    palette <- getPalatte(length(unique(df.union$Clonotype)))
  }
  g <- ggplot(df.union, aes(x=round, y=Count)) +
    geom_bar(stat='identity', aes(fill=Clonotype)) +
    theme(legend.position="bottom", legend.box = "horizontal", legend.title=element_blank(), legend.text=element_text(size=7)) +
    labs(title=paste("Top", top, "clonotype across each sample"),
         subtitle=paste0("Colour coded clonotypes, distribution of each clonotype is relative to top ", top, " - not overall."),
         x="Sample",
         y="Distribution")  +
    scale_fill_manual(values=palette)
  return (g)
}

#library(ggplot2)
#library(RColorBrewer)
#fs <- list.files(pattern = "PCR[123].*_cdr3_clonotypes_all_over.csv.gz", recursive = TRUE, full.names = TRUE)
#p <- topNDist(lapply(fs, read.csv, stringsAsFactors = FALSE), c("PCR1", "PCR2", "PCR3"))
#fs <- list.files(pattern = "PCR[12345].*_cdr3_clonotypes_all_over.csv.gz", recursive = TRUE, full.names = TRUE)
#p <- topNDist(lapply(fs, read.csv, stringsAsFactors = FALSE), c("PCR1", "PCR2", "PCR3", "PCR4", "PCR5"))
#plot(p)
