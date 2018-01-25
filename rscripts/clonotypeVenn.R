library(VennDiagram)

vennIntersection <- function(dataframes, sampleNames, outFile, top) {
  # Plots a venn diagram for (non-)intersecting clonotypes
  # Args:
  #     dataframes: A list() type. dataframes of sample. Accepts from 2 - 5 samples. (incl) Will not plot for anything else
  #     sampleNames: A vector type. Vector of strings, each representing sample name
  #     outFile: A string type. fully qualified path with output file name (assumes PNG for now)
  #     top: Top N clonotypes to plot. If unspecified, defaults to ALL clones available in dataframes
  # Returns: Nothing. Plots and saves venndiagram
  nsample <- length(dataframes)

  stopifnot(nsample == length(sampleNames))
  if (nsample >= 2 && nsample <= 5) {
      # output
      png(file = outFile, width = 8, height = 7, units = "in", res = 300)
      
      # cleanups
      colNames <- c("Clonotype", "Count")
      for (i in 1:nsample) {
        df <- dataframes[[i]]
        if (!missing(top)) {
          df <- head(df[colNames], top)
        } else {
          df <- df[colNames]
        }
        dataframes[[i]] <- df
      }
      
      # merge all dataframes
      df.union <- merge(dataframes[[1]], dataframes[[2]], by = "Clonotype", all.y = TRUE, all.x = TRUE)
      colnames(df.union) <- c("Clonotype", sampleNames[1], sampleNames[2])
      df.union[is.na(df.union)] <- 0
      # check if there's more
      if (nsample > 2) {
        for (i in 3:nsample) {
          oldColNames <- colnames(df.union)
          df.union <- merge(df.union, dataframes[[i]], by = "Clonotype", all.y = TRUE, all.x = TRUE)
          colnames(df.union) <- c(oldColNames, sampleNames[i])
          df.union[is.na(df.union)] <- 0
        }
      }
      
      # plot!
      area1 <- sum(df.union[sampleNames[1]] > 0)
      area2 <- sum(df.union[sampleNames[2]] > 0)
      area12 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0)
      if (nsample == 2) {
        draw.pairwise.venn(area1, area2, area12, category = c(sampleNames[1], sampleNames[2]),
                           lty = "blank",
                           col = "transparent",
                           cat.fontface = "bold",
                           cex = 1.7,
                           fill = 2:3,
                           alpha = 0.5,
                           filename = NULL)
      } else {
        area3 <- sum(df.union[sampleNames[3]] > 0)
        area13 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[3]] > 0)
        area23 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0)
        area123 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0)
        if (nsample == 3) {
          draw.triple.venn(area1 = area1, area2 = area2, area3 = area3, n12 = area12, n13 = area13, n23 = area23, n123 = area123,
                           category = c(sampleNames[1], sampleNames[2], sampleNames[3]),
                           lty = "blank",
                           col = "transparent", cat.fontface = "bold", cex = 1.7,
                           fill = 2:4, alpha = 0.5, filename = NULL)
        } else {
            area4 <- sum(df.union[sampleNames[4]] > 0)
            area14 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[4]] > 0)
            area24 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[4]] > 0)
            area34 <- sum(df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0)
            area124 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0 & df.union[sampleNames[4]] > 0)
            area134 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0)
            area234 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0)
            area1234 <- sum(df.union[sampleNames[1]] > 0 &
                            df.union[sampleNames[2]] > 0 &
                            df.union[sampleNames[3]] > 0 &
                            df.union[sampleNames[4]] > 0)
            if (nsample == 4) {
                draw.quad.venn(area1 = area1, area2 = area2, area3 = area3, area4 = area4,
                               n12 = area12, n13 = area13, n14 = area14, n23 = area23, n24 = area24, n34 = area34,
                               n123 = area123, n124 = area124, n134 = area134, n234 = area234, n1234 = area1234,
                               category = c(sampleNames[1], sampleNames[2], sampleNames[3], sampleNames[4]),
                               lty = "blank",
                               col = "transparent", cat.fontface = "bold", cex = 1.7,
                               fill = 2:5, alpha = 0.5, filename = NULL
                )
            } else {
                # quuntuple plot
                area5 <- sum(df.union[sampleNames[5]] > 0)
                area15 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[5]] > 0)
                area25 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[5]] > 0)
                area35 <- sum(df.union[sampleNames[3]] > 0 & df.union[sampleNames[5]] > 0)
                area45 <- sum(df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                area125 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[2]] > 0 & df.union[sampleNames[5]] > 0)
                area135 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[5]] > 0)
                area145 <- sum(df.union[sampleNames[1]] > 0 & df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                area235 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[3]] > 0 & df.union[sampleNames[5]] > 0)
                area245 <- sum(df.union[sampleNames[2]] > 0 & df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                area345 <- sum(df.union[sampleNames[3]] > 0 & df.union[sampleNames[4]] > 0 & df.union[sampleNames[5]] > 0)
                area1235 <- sum(df.union[sampleNames[1]] > 0 &
                                df.union[sampleNames[2]] > 0 &
                                df.union[sampleNames[3]] > 0 &
                                df.union[sampleNames[5]] > 0)
                area1245 <- sum(df.union[sampleNames[1]] > 0 &
                                df.union[sampleNames[2]] > 0 &
                                df.union[sampleNames[4]] > 0 &
                                df.union[sampleNames[5]] > 0)
                area1345 <- sum(df.union[sampleNames[1]] > 0 &
                                df.union[sampleNames[3]] > 0 &
                                df.union[sampleNames[4]] > 0 &
                                df.union[sampleNames[5]] > 0)
                area2345 <- sum(df.union[sampleNames[2]] > 0 &
                                df.union[sampleNames[3]] > 0 &
                                df.union[sampleNames[4]] > 0 &
                                df.union[sampleNames[5]] > 0)
                area12345 <- sum(df.union[sampleNames[1]] > 0 &
                                 df.union[sampleNames[2]] > 0 &
                                 df.union[sampleNames[3]] > 0 &
                                 df.union[sampleNames[4]] > 0 &
                                 df.union[sampleNames[5]] > 0)
                draw.quintuple.venn(
                               area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5,
                               n12 = area12, n13 = area13, n14 = area14, n15 = area15, n23 = area23, n24 = area24, n25 = area25, n34 = area34,
                               n35 = area35, n45 = area45, n123 = area123, n124 = area124, n125 = area125, n134 = area134, 
                               n135 = area135, n145 = area145, n234 = area234, n235 = area235, n245 = area245, n345 = area345, n1234 = area1234,
                               n1235 = area1235, n1245 = area1245, n1345 = area1345, n2345 = area2345, n12345 = area12345,
                               category = c(sampleNames[1], sampleNames[2], sampleNames[3], sampleNames[4], sampleNames[5]),
                               lty = "blank",
                               col = "transparent", cat.fontface = "bold", cex = 1.7,
                               fill = 2:6, alpha = 0.5, filename = NULL
                )
            }
            
        }
      }
      dev.off()
  }
}

#fs <- list.files(pattern = "PCR[123].*_cdr3_clonotypes_all_over.csv.gz", recursive = TRUE, full.names = TRUE)
#vennIntersection(lapply(fs, read.csv, stringsAsFactors = FALSE), c("PCR1", "PCR2", "PCR3"), "./venn.png")
