# single plot uses blue
BLUEHEX <- "#56B4E9"



checkVert <- function(filename) {
  f <- file(filename, "r")
  res <- readLines(f, n = 1) == "vert"
  close(f)
  return (res)
}

listFilesInOrder <- function(path, pattern, expectedRet = 1) {
  # Returns files in order of path. I.e overrides list.files' behaviour of sorted files
  # This is crucial because we are assuming sampleNames is in 1-1 correspondance with dataframes.
  # This can be achieved by manually iterating over all sample path in the provided vector of path
  # and appending to a vector. Sometimes, (like abundance) list.files will return more that one
  # matching file (gene, family, variant), so expectedRet is there to ensure that we aren't doing
  # something silly.
  # Returns: ordered vector of files (according to provieded path's ordering)
  orderedFiles <- c()
  for (p in path) {
    retval <- list.files(path = p, pattern = pattern, full.names = TRUE, recursive = TRUE)
    stopifnot(length(retval) == expectedRet)
    orderedFiles <- c(orderedFiles,  retval)
  }
  return (orderedFiles)
}