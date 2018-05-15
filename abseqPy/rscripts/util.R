# single plot uses blue
BLUEHEX <- "#56B4E9"

# plot sizes
V_WIDTH <- 8
V_HEIGHT <- 5
H_WIDTH <- 5
H_HEIGHT <- 8
VENN_WIDTH <- 8
VENN_HEIGHT <- 8
V_WIDTH_L <- 12
V_HEIGHT_L <- 7.5



checkVert <- function(filename) {
  f <- file(filename, "r")
  res <- grepl("vert", readLines(f, n = 1), fixed = TRUE)
  close(f)
  return (res)
}

getTotal <- function(filename) {
  f <- file(filename, "r")
  res <- unlist(strsplit(readLines(f, n = 1), "="))[2]
  close(f)
  return (res)
}

listFilesInOrder <- function(path, pattern, expectedRet = c(1)) {
  # Returns files in order of path. I.e overrides list.files' behaviour of sorted files
  # This is crucial because we are assuming sampleNames is in 1-1 correspondance with dataframes.
  # This can be achieved by manually iterating over all sample path in the provided vector of path
  # and appending to a vector. Sometimes, (like abundance) list.files will return more that one
  # matching file (gene, family, variant), so expectedRet is there to ensure that we aren't doing
  # something silly. expectRet is a vector because sometimes there might be more than one possible configuration.
  # EG: abundance plot may or maynot have D gene when analyzing heavy/light chains
  # Returns: ordered vector of files (according to provided path's ordering)
  orderedFiles <- c()
  for (p in path) {
    retval <- list.files(path = p, pattern = pattern, full.names = TRUE, recursive = TRUE)
    if (length(retval) == 0) {
        return (c())
    }
    stopifnot(length(retval) %in% expectedRet)
    orderedFiles <- c(orderedFiles,  retval)
  }
  return (orderedFiles)
}

inferAnalyzed <- function(sampleDirectory) {
  return (list.files(sampleDirectory))
}

capitalize <- function(str) {
  firstLetter <- substr(str, 1, 1)
  rest <- substr(str, 2, nchar(str))
  return (paste0(toupper(firstLetter), rest))
}