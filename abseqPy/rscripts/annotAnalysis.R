annotAnalysis <- function(annotDirectories, annotOut, sampleNames, mashedNames) {
  # with outliers
  searchFiles <- listFilesInOrder(path = annotDirectories, pattern = ".*_all_clones_len_dist\\.csv(\\.gz)?$")
  if (length(searchFiles) > 0) {
      g <- plotSpectratype(
        lapply(searchFiles, read.csv),
        sampleNames,
        title = "Sequence lengths",  # plotSpectratype names the sample(s) for us in the title; don't have to do it here
        xlabel = "Sequence Length (bp)",
        ylabel = "Proportion"
      )
      ggsave(paste0(annotOut, mashedNames, "_all_clones_len_dist.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  }
  # without outliers
  searchFiles <- listFilesInOrder(path = annotDirectories, pattern = ".*_all_clones_len_dist_no_outliers\\.csv(\\.gz)?$")
  if (length(searchFiles) > 0) {
      g <- plotSpectratype(
        lapply(searchFiles, read.csv),
        sampleNames,
        title = "Sequence lengths",  # plotSpectratype names the sample(s) for us in the title; don't have to do it here
        xlabel = "Sequence Length (bp)",
        ylabel = "Proportion"
      )
      ggsave(paste0(annotOut, mashedNames, "_all_clones_len_dist_no_outliers.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
  }
}