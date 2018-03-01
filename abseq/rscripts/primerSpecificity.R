allPrimerNames <- function(primerFile) {
  primers <- c()
  fp = file(primerFile, open="r")
  while (TRUE) {
    line = readLines(fp, n = 1)
    if (!length(line)) {
      break
    }
    if (startsWith(line, ">")) {
      primers <- c(primers, gsub("^>\\s*", "", line))
    }
  }
  close(fp)
  return (primers)
}

canonicalizeTitle <- function(str) {
  if (str == "outframe") {
    return ("Out-of-frame")
  } else if (str == "indel_pos") {
    return ("Abundance of Indel Positions")
  } else if (str == "indelled") {
    return ("Indelled abundance") 
  } else {
    return (capitalize(str))
  }
}

primerAnalysis <- function(primer5File, primer3File, primerDirectories, primerOut, combinedNames, mashedNames) {
  if (primer5File == "None") {
    primer5 <- c()
  } else {
    primer5 <- allPrimerNames(primer5File)
  }
  if (primer3File == "None") {
    primer3 <- c()
  } else {
    primer3 <- allPrimerNames(primer3File)
  }
  
  allPrimers <- list(primer5, primer3)
  category <- c("all", "productive", "outframe")
  analysisType <- c("stopcodon", "integrity", "indelled", "indel_pos")
  
  for (c in 1:length(category)) {
    
    for (i in 1:length(allPrimers)) {
      # what end this is, 3? 5?
      primerNames <- allPrimers[[i]]
      if (length(primerNames)) {
        if (i == 1) {
          pend <- '5'
        } else {
          pend <- '3'
        }
        
        ####################################################################################
        #                           individual primer IGV abundance                        #
        ####################################################################################
        for (j in 1:length(primerNames)) {
          files <- listFilesInOrder(path = primerDirectories, pattern = paste0(".*_", category[c], "_", pend, "end_", primerNames[j],"_igv_dist\\.csv(\\.gz)?$"))
          # there is a slight chance that the user provided primer has no hit, hence, not dist files (AbSEq (obviously) doesn't plot something empty)
          if (length(files)) {
            subtitle <- paste("Total is", paste(lapply(files, function(x) { as.integer(getTotal(x)) }), collapse = ', '))
            vertical<- checkVert(files[[1]])
            primPlot <- plotDist(
              lapply(files, read.csv, skip = 1),
              sampleNames,
              paste(paste0("IGV Abundance of ", primerNames[j]," (", canonicalizeTitle(category[c]), ") in "), combinedNames),
              vertical,
              subs = subtitle
            )
            if (vertical) {
              ggsave(paste0(primerOut, mashedNames, paste0("_", category[c], "_", pend, "end_", primerNames[j], "_igv_dist.png")), plot = primPlot, width = V_WIDTH, height = V_HEIGHT)
            } else {
              ggsave(paste0(primerOut, mashedNames, paste0("_", category[c], "_", pend, "end_", primerNames[j], "_igv_dist.png")), plot = primPlot, width = H_WIDTH, height = H_HEIGHT)
            }
          }  # endof if length(files)
        }   # end of individual primer IGV abundance
        
        
        ####################################################################################
        #      stopcodon, integrity, indelled and indel_pos analysis on all categories     #
        ####################################################################################
        for (j in 1:length(analysisType)) {
          files <- listFilesInOrder(path = primerDirectories, pattern = paste0(".*_", category[c], "_", pend, "end_", analysisType[j], "_dist\\.csv(\\.gz)?$"))
          # there is a slight chance that the user provided primer has no hit, hence, not dist files (AbSEq (obviously) doesn't plot something empty)
          if (length(files)) {
            subtitle <- paste("Total is", paste(lapply(files, function(x) { as.integer(getTotal(x)) }), collapse = ', '))
            vertical <- checkVert(files[[1]])
            primPlot <- plotDist(
              lapply(files, read.csv, skip = 1),
              sampleNames,
              paste(paste0(canonicalizeTitle(analysisType[j]), " of ", pend, "'-end Primer Sequence (", canonicalizeTitle(category[c]), ") in"), combinedNames),
              vertical,
              subs = subtitle
            )
            if (vertical) {
              ggsave(paste0(primerOut, mashedNames, paste0("_", category[c], "_", pend, "end_", analysisType[j], "_dist.png")), plot = primPlot, width = V_WIDTH, height = V_HEIGHT)
            } else {
              ggsave(paste0(primerOut, mashedNames, paste0("_", category[c], "_", pend, "end_", analysisType[j], "_dist.png")), plot = primPlot, width = H_WIDTH, height = H_HEIGHT)
            }
          }  # endof if length(files)
        } # endof "type" analysis
      } # endof if length(primerNames)         
    } # end of for (i in allPrimers)
  } # end of category
}
