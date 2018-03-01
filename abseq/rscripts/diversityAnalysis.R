diversityAnalysis <- function(diversityDirectories, diversityOut, sampleNames, mashedNames) {
    #               CLONOTYPE PLOTS                  #
    ##################################################
    cdr3ClonesFile <- listFilesInOrder(path = diversityDirectories, pattern = ".*_cdr3_clonotypes_.*_over\\.csv(\\.gz)?$")
    if (length(cdr3ClonesFile) > 0 && length(sampleNames) > 1) {
      # plot scatter plot (CDR3 clonotypes)
      scatterClones(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
                    sampleNames,
                    diversityOut, "CDR3")
      
      # plot venn diagram (clonotypes)
      vennIntersection(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE), 
                       sampleNames,
                       paste0(diversityOut, mashedNames, "_cdr3_clonotypeIntersection.png"))
      
      # plot top N distribution (clonotypes)
      g <- topNDist(lapply(cdr3ClonesFile, read.csv, stringsAsFactors = FALSE),
                    sampleNames)
      ggsave(paste0(diversityOut, mashedNames, "_top10Clonotypes.png"), plot = g, width = V_WIDTH_L, height = V_HEIGHT_L)
    }
    
    #               FR/CDR    PLOTS                  #
    ##################################################
    # plot duplication, rarefaction, recapture
    for (reg in c("cdr", "cdr_v", "fr")) {
      if (reg == "cdr")
        includedRegions <- c("CDR1", "CDR2", "CDR3")
      else if (reg == "cdr_v")
        includedRegions <- c("CDR3", "V")
      else
        includedRegions <- c("FR1", "FR2", "FR3")
      searchFiles <- listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_", reg,"_duplication\\.csv(\\.gz)?$"))
      if (length(searchFiles) > 0) {
          # duplication
          g <- plotDuplication(searchFiles,
                               sampleNames,
                               includedRegions)
          ggsave(paste0(diversityOut, mashedNames, "_", reg, "_duplication.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
      }
      
      searchFiles <- listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_", reg, "_rarefaction\\.csv(\\.gz)?$"))
      if (length(searchFiles) > 0) {
          # plot rarefaction
          g <- plotRarefaction(searchFiles,
                               sampleNames,
                               includedRegions)
          ggsave(paste0(diversityOut, mashedNames, "_", reg, "_rarefaction.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
      }
      searchFiles <- listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_", reg, "_recapture\\.csv(\\.gz)?$"))
      if (length(searchFiles) > 0) {
          # plot recapture
          g <- plotRecapture(searchFiles,
                             sampleNames,
                             includedRegions)
          ggsave(paste0(diversityOut, mashedNames, "_", reg, "_recapture.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
      }
    }
    
    # MINISECTION:
    # ## SPECTRATYPES ###
    specOut <- paste0(diversityOut, "spectratypes/")
    dir.create(specOut)
    # CDR 1 - 3
    for (i in 1:3) {
      specFiles <- listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_cdr", i, "_spectratype\\.csv(\\.gz)?$"))
      if (length(specFiles) > 0) {
          g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                               sampleNames,
                               paste0("CDR", i))
          ggsave(paste0(specOut, mashedNames, "_cdr", i, "_spectratype.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
      }
    }
    
    # special case, no outliers plot for CDR3 only
    specFiles <- listFilesInOrder(path = diversityDirectories, pattern = ".*_cdr3_spectratype_no_outliers\\.csv(\\.gz)?$")
    if (length(specFiles) > 0) {
        g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                             sampleNames,
                             "CDR3")
        ggsave(paste0(specOut, mashedNames, "_cdr3_spectratype_no_outliers.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
    }
    
    # FR 1 - 4
    for (i in 1:4) {
      specFiles <- listFilesInOrder(path = diversityDirectories, pattern = paste0(".*_fr", i, "_spectratype\\.csv(\\.gz)?$"))
      if (length(specFiles) > 0) {
          g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                               sampleNames,
                               paste0("FR", i))
          ggsave(paste0(specOut, mashedNames, "_fr", i, "_spectratype.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
      }
    }
    
    # entire V-domain
    specFiles <- listFilesInOrder(path = diversityDirectories, pattern = ".*_v_spectratype\\.csv(\\.gz)?$")
    if (length(specFiles) > 0) {
        g <- plotSpectratype(lapply(specFiles, read.csv, stringsAsFactors = FALSE),
                             sampleNames,
                             "V domain")
        ggsave(paste0(specOut, mashedNames, "_v_spectratype.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
    }
    
    # we can plot region analysis if there's only one sample
    #if (length(sampleNames) == 1) {
    #  # default = top 15
    #  g <- regionAnalysis(diversityOut, sampleNames[1])
    #  ggsave(paste0(diversityOut, mashedNames, "_region_analysis.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
    #}
}