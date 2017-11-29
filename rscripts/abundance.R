library(ggplot2)

abundancePlot <- function(fs, sampleNames, outputDir) {
  # Plots 6 abundance plots
  # Args:
  #     fs: A List() type. List of files with VDJ distributions (under _ig[VDJ]_dist_[family|variant|gene]_level.csv)
  #     sampleNames: Vector type. Vector of strings, each representing sample name (1-1 with fs)
  #     outputDir: Output directory for png files.
  for (expression in c("family", "gene")) {
    for (gene in c('v', 'd', 'j')) {
      # correction, J has no "gene" but rather variant
      if (gene == 'j' && expression == 'gene') {
        expression <- 'variant'
      }
      reg <- paste(".*_ig", gene, "_dist_", expression, "_level\\.csv(\\.gz)?$", sep = "")
      selectedFiles <- fs[grepl(reg, fs)]
      vert <- checkVert(selectedFiles[1])
      mashedName <- paste(sampleNames, collapse = ", ")
      dataframes <- lapply(selectedFiles, read.csv, stringsAsFactors=FALSE, skip = 1)
      p <- plotDist(dataframes, sampleNames, paste("IG", toupper(gene), " abundance in ", mashedName, sep = ""), vert)
      ggsave(paste(outputDir, paste(sampleNames, collapse = "_"), "_ig", gene, "_dist_", expression, "_level.png", sep=""), plot = p)
    }
  }
}

# https://stackoverflow.com/questions/17499013/how-do-i-make-a-list-of-data-frames
# test driver:
# fs can be obtained as such fs <- list.files(pattern = "\\.csv$")
#fs <- list.files(pattern = "PCR[123].*ig[vdj]_dist_[family|gene|variant].*\\.csv$", full.names = TRUE, recursive = TRUE)
#abundancePlot(fs, c("PCR1", "PCR2", "PCR3"), "/Users/harry/")
