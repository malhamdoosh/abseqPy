upstreamAnalysis <- function(upstreamDirectories, upstreamOut, expectedLength, upstreamLengthRange, sampleNames, combinedNames, mashedNames, secsig) {
  if (secsig) {
    fname <- "secsig"
  } else {
    fname <- "5utr"
  }
    
  if (!is.na(expectedLength)) {
    level <- c("gene", "family")
    for (lvl in level) {
      statuses <- c("valid", "faulty", "no_atg")
      for (status in statuses) {
        requestedFiles <- listFilesInOrder(path = upstreamDirectories,
                                           pattern = paste0(".*_", fname, "_", upstreamLengthRange, "_", status,"_", lvl, '\\.csv(\\.gz)?$')
                                           )
        if (length(requestedFiles) > 0) {
          if (status == "valid") {
            if (secsig) {
              title <- "Valid Secretion Signals"
            } else {
              title <- "Valid 5'-UTRs"
            }
          } else if (status == 'faulty') {
            title <- "Faulty Translations"
          } else {
            title <- "Upstream sequences without start codon"
          }
          subtitle <- paste("Total is ", paste(lapply(requestedFiles, function(x) { as.integer(getTotal(x)) }), collapse = ", "))
          plotVert <- checkVert(requestedFiles[[1]])
          g <- plotDist(lapply(requestedFiles, read.csv, skip = 1), sampleNames, title, vert = plotVert, subs = subtitle)
          if (plotVert) {
            width <- V_WIDTH;
            height <- V_HEIGHT;
          } else {
            width <- H_WIDTH;
            height <- H_WIDTH;
          }
          ggsave(paste0(upstreamOut,
                        mashedNames,
                        paste0("_", fname, "_", upstreamLengthRange, "_", status, "_", lvl, ".png")),
                 plot = g,
                 width = width,
                 height = height)
        }
      }
    }
  }
  
}

upstreamDist <- function(upstreamDirectories, upstreamOut, exepctedLength, upstreamLengthRange, sampleNames, combinedNames, mashedNames, secsig) {
  if (secsig) {
    fname <- "secsig"
  } else {
    fname <- "5utr"
  }
    
  # full lengthed upstream sequences and upstream seqs that are shorter than the expected length
  if (!is.na(expectedLength)) {
    lengths <- c("", "_short")
  } else {
    # if expected length is NA => Expected length is Inf
    lengths <- c("")
  }
  
  for (len in lengths) {
    seqLengthFiles <- listFilesInOrder(path = upstreamDirectories, pattern = paste0(".*_", fname, "_", upstreamLengthRange, "_dist", len, "\\.csv(\\.gz)?$"))
    if (length(seqLengthFiles) > 0) {
      g <- plotSpectratype(lapply(seqLengthFiles, read.csv),
                           sampleNames,
                           title = "Sequence lengths",
                           xlabel = "Sequence Length(bp)",
                           ylabel = "Distribution")
      ggsave(paste0(upstreamOut, mashedNames, paste0("_", fname, "_", upstreamLengthRange, "_dist", len, ".png")), plot = g, width = V_WIDTH, height = V_HEIGHT)
    }
    
    # class level
    seqClassLengthFiles <- listFilesInOrder(path = upstreamDirectories, pattern = paste0(".*_", fname, "_", upstreamLengthRange,"_dist", len,"_class", "\\.csv(\\.gz)?$"))
    if (length(seqClassLengthFiles) > 0) {
      subtitle <- paste("Total is ", paste(lapply(seqClassLengthFiles, function(x) { as.integer(getTotal(x)) }), collapse = ", "))
      plotVert <- checkVert(seqClassLengthFiles[[1]])
      g <- plotDist(lapply(seqClassLengthFiles, read.csv, skip = 1), 
                    sampleNames,
                    paste("IGV Abundance in Sample", combinedNames),
                    plotVert,
                    subs = subtitle)
      if (plotVert) {
        ggsave(paste0(upstreamOut, mashedNames, paste0("_", fname, "_", upstreamLengthRange, "_dist", len, "_class",".png")), plot = g, width = V_WIDTH, height = V_HEIGHT)
      } else {
        ggsave(paste0(upstreamOut, mashedNames, paste0("_", fname, "_", upstreamLengthRange, "_dist", len, "_class",".png")), plot = g, width = H_WIDTH, height = H_HEIGHT)
      }
    }
    
    # box plot for class level
    seqClassLengthBoxFiles <- listFilesInOrder(path = upstreamDirectories, pattern = paste0(".*_", fname, "_", upstreamLengthRange, "_dist", len, "_class", "_box\\.csv(\\.gz)?$"))
    if (length(seqClassLengthBoxFiles) > 0) {
      g <- boxPlot(lapply(seqClassLengthBoxFiles, read.csv), sampleNames, paste("Sequence Lengths in", combinedNames))
      ggsave(paste0(upstreamOut, mashedNames, paste0("_", fname, "_", upstreamLengthRange, "_dist", len, "_class", "_box.png")), plot = g, width = V_WIDTH, height = V_HEIGHT)
    }
  }
}