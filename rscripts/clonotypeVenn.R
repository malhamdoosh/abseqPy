library(VennDiagram)

vennIntersection <- function(dataframes, sampleNames, outFile, top) {
  # Plots a venn diagram for (non-)intersecting clonotypes
  # Args:
  #     dataframes: A list() type. dataframes of sample. As of 20/11/2017, only supports nsample = 2 or 3
  #     sampleNames: A vector type. Vector of strings, each representing sample name
  #     outFile: A string type. fully qualified path with output file name (assumes PNG for now)
  #     top: Top N clonotypes to plot. If unspecified, defaults to ALL clones available in dataframes
  # Returns: Nothing. Plots and saves venndiagram
  nsample <- length(dataframes)
  stopifnot(nsample > 1 && nsample <= 3) # TODO: extend this to 6
  stopifnot(nsample == length(sampleNames))
  
  # output
  png(file = outFile)
  
  # cleanups
  colNames <- c("Clonotype", "Count")
  for (i in 1:nsample) {
    df <- dataframes[[i]]
    if (!missing(top)) {
      df <- head(df[colNames], top)
    } else {
      df <- df[colNames]
    }
    dataframes[[i]] <- df
  }
  
  # merge all dataframes
  df.union <- merge(dataframes[[1]], dataframes[[2]], by = "Clonotype", all.y = TRUE, all.x = TRUE)
  colnames(df.union) <- c("Clonotype", sampleNames[1], sampleNames[2])
  df.union[is.na(df.union)] <- 0
  # check if there's more
  if (nsample > 2) {
    for (i in 3:nsample) {
      oldColNames <- colnames(df.union)
      df.union <- merge(df.union, dataframes[[i]], by = "Clonotype", all.y = TRUE, all.x = TRUE)
      colnames(df.union) <- c(oldColNames, sampleNames[i])
      df.union[is.na(df.union)] <- 0
    }
  }
  
  # plot!
  area1 <- sum(df.union[sampleNames[1]] > 0)
  area2 <- sum(df.union[sampleNames[2]] > 0)
  area12 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0)
  if (nsample == 2) {
    draw.pairwise.venn(area1, area2, area12, category = c(sampleNames[1], sampleNames[2]),
                       lty = "blank",
                       col = "transparent",
                       cat.fontface = "bold",
                       cex = 1.7,
                       fill = 2:3,
                       alpha = 0.5,
                       filename = NULL)
  } else {
    area3 <- sum(df.union[sampleNames[3]] > 0)
    area13 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[3]] > 0)
    area23 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0)
    area123 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0)
    if (nsample == 3) {
      draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, n12 = area12, n13 = area13, n23 = area23, n123 = area123,
                       category = c(sampleNames[1], sampleNames[2], sampleNames[3]),
                       lty = "blank",
                       col = "transparent", cat.fontface = "bold", cex = 1.7,
                       fill = 2:4, alpha = 0.5, filename = NULL)
    } else {
      # not supported yet! (TODO)
    }
  }
  dev.off()
}

#fs <- list.files(pattern = "PCR[123].*_cdr3_clonotypes_all_over.csv.gz", recursive = TRUE, full.names = TRUE)
#vennIntersection(lapply(fs, read.csv, stringsAsFactors = FALSE), c("PCR1", "PCR2", "PCR3"), "/Users/harry/venn.png")