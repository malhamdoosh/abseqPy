productivityAnalysis <- function(productivityDirectories, prodOut, sampleNames, combinedNames, mashedNames) {
    # main productivity file
    ########################
    prodFiles <- listFilesInOrder(path = productivityDirectories, pattern = ".*_productivity\\.csv(\\.gz)?$")
    if (length(prodFiles) > 0) {
        g <- productivityPlot(lapply(prodFiles, read.csv, stringsAsFactors = FALSE), sampleNames)
        ggsave(paste0(prodOut, mashedNames, "_productivity.png"), plot = g, width = V_WIDTH, height = V_HEIGHT)
    }
    
    # sub-productivity files
    ########################
    regions <- c("cdr1", "cdr2", "cdr3", "fr1", "fr2", "fr3", "igv", "igd", "igj")
    
    # gaps_dist plots only
    gapPlots <- prodDistPlot(productivityDirectories,
                              sampleNames,
                              "Gaps in",
                              "_gaps_dist\\.csv(\\.gz)?$",
                              unlist(lapply(regions, function(region) { paste0(prodOut, mashedNames, "_", region, "_gaps_dist.png") })),
                              regions)
    # gaps_out_of_frame plots only (no igv, igd, ihj plots for this)
    subregions <- head(regions, n = 6)
    gapOutFramePlots <- prodDistPlot(productivityDirectories,
                                     sampleNames,
                                     "Gaps in",
                                     "_gaps_dist_out_of_frame\\.csv(\\.gz)?$",
                                     unlist(lapply(subregions, function(region) {paste0(prodOut, mashedNames, "_", region, "_gaps_dist_out_of_frame.png") })),
                                     subregions)
    # mismatch dist only
    mismatchPlots <- prodDistPlot(productivityDirectories,
                                  sampleNames,
                                  "Mismatches in",
                                  "_mismatches_dist\\.csv(\\.gz)?$",
                                  unlist(lapply(regions, function(region) { paste0(prodOut, mashedNames, "_", region, "_mismatches_dist.png")})),
                                  regions)
    # stop codon dist plot
    stopCodonPlot <- prodDistPlot(productivityDirectories, sampleNames,
                              "Stop codon in In-frame Clones",
                              "_stopcodon_dist_in_frame\\.csv(\\.gz)?$",
                              c(paste0(prodOut, mashedNames, "_stopcodon_dist_in_frame.png")),
                              c("")) # no regions
    # vjframe plot
    vjframePlot <- prodDistPlot(productivityDirectories, sampleNames,
                                "V-D-J Rearrangement",
                                "_vjframe_dist\\.csv(\\.gz)?$",
                                c(paste0(prodOut, mashedNames, "_vjframe_dist.png")),
                                c("")) # no regions
    # 3 special cases for IGV region
    # igv - inframe_unproductive, out of frame, productive
    igvInframeUnProductive <- prodDistPlot(productivityDirectories, sampleNames,
                                        "Abundance of In-frame Unproductive Clones in",
                                        "_dist_inframe_unproductive\\.csv(\\.gz)?$",
                                        c(paste0(prodOut, mashedNames, "_igv_dist_inframe_unproductive.png")),
                                        c("igv"))
    
    igvOutOfFrame <- prodDistPlot(productivityDirectories, sampleNames,
                                  "Abundance of Out-Of-Frame Clones in",
                                  "_dist_out_of_frame\\.csv(\\.gz)?$",
                                  c(paste0(prodOut, mashedNames, "_igv_dist_out_of_frame.png")),
                                  c("igv"))
    
    igvProductivity <- prodDistPlot(productivityDirectories, sampleNames,
                                    "Abundance of Productive Clones in",
                                    "_dist_productive\\.csv(\\.gz)?$",
                                    c(paste0(prodOut, mashedNames, "_igv_dist_productive.png")),
                                    c("igv"))
    
    # stop codon in FR/CDR region proportion plot
    #############################################
    frameStatus <- c("inframe", "outframe")
    for (framestat in frameStatus) {
      stopCodonRegionFiles <- listFilesInOrder(path = productivityDirectories, pattern = paste0(".*_stopcodon_region_", framestat, "\\.csv(\\.gz)?$"))
      if (length(stopCodonRegionFiles) > 0) {
        vert = checkVert(stopCodonRegionFiles[[1]])
        if (framestat == 'inframe') {
          titlestatus <- "In-frame"
        } else {
          titlestatus <- "Out-of-frame"
        }
        subtitle <- paste("Total is", paste(lapply(stopCodonRegionFiles, function(x) {as.integer(getTotal(x)) }), collapse = ", "))
        stopcodonRegion <- plotDist(
          lapply(stopCodonRegionFiles, read.csv, skip = 1),
          sampleNames,
          paste("Stop codon in FRs and CDRs of", titlestatus, "sequences in", combinedNames),
          vert,
          subs = subtitle,
          sortDist = FALSE
        )
        if (vert) {
          ggsave(paste0(prodOut, mashedNames, "_stopcodon_region_", framestat, ".png"), plot = stopcodonRegion, width = V_WIDTH, height = V_HEIGHT)
        } else {
          ggsave(paste0(prodOut, mashedNames, "_stopcodon_region_", framestat, ".png"), plot = stopcodonRegion, width = H_WIDTH, height = H_HEIGHT)
        }
      }
    }
}