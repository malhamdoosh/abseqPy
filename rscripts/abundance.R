library(ggplot2)

abundancePlot <- function(fs, sampleNames, outputDir) {
  # Plots 6 abundance plots (V, D) => gene, family, (J) => variant, family
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
      reg <- paste0(".*_ig", gene, "_dist_", expression, "_level\\.csv(\\.gz)?$")
      selectedFiles <- fs[grepl(reg, fs)]
      if (length(selectedFiles) > 0) {
          vert <- checkVert(selectedFiles[1])
          if (vert) {
            width <- V_WIDTH
            height <- V_HEIGHT
          } else {
            width <- H_WIDTH
            height <- H_HEIGHT
          }
          mashedName <- paste(sampleNames, collapse = ", ")
          dataframes <- lapply(selectedFiles, read.csv, stringsAsFactors=FALSE, skip = 1)
          subs <- "Total is"
          frames <- length(dataframes)
          subs <- paste(subs, paste(lapply(selectedFiles, function(x) { round(getTotal(x)) }), collapse = ", "))
          p <- plotDist(dataframes, sampleNames, paste0("IG", toupper(gene), " abundance in ", mashedName), vert, subs = subs)
          ggsave(paste0(outputDir, paste(sampleNames, collapse = "_"), "_ig", gene, "_dist_", expression, "_level.png"), plot = p, width = width, height = height)
      }
    }
  }
}

# https://stackoverflow.com/questions/17499013/how-do-i-make-a-list-of-data-frames
# test driver:
# fs can be obtained as such fs <- list.files(pattern = "\\.csv$")
#fs <- list.files(pattern = "PCR[123].*ig[vdj]_dist_[family|gene|variant].*\\.csv$", full.names = TRUE, recursive = TRUE)
#abundancePlot(fs, c("PCR1", "PCR2", "PCR3"), "/Users/harry/")
