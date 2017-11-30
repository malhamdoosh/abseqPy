# single plot uses blue
BLUEHEX <- "#56B4E9"



checkVert <- function(filename) {
  f <- file(filename, "r")
  res <- readLines(f, n = 1) == "vert"
  close(f)
  return (res)
}

