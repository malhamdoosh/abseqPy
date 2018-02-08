library(circlize)
library(RColorBrewer)

plotCirclize <- function(sampleName, path) {
  # Plots a circlize chord Diagram for the V-J association of a given sample
  # 
  # Args:
  #    sampleName: A string type. Name of sample
  #    path: A string type. Path to vjassoc.csv directory
  #
  # Returns: Nothing. Saves a png named <sampleName>_vjassoc.csv in the same directory as path
  
  # open and read csv
  filename <- paste0(path, sampleName, "_vjassoc.csv")
  print(filename)
  if (file.exists(filename)) {
      df <- read.csv(filename)
      
      # output file
      png(gsub(".csv", ".png", filename))
      
      # circos theme setup
      if (length(unique(df[[1]]))-1 < 10 && length(unique(df[[2]]))-1 < 10)  {
          circos.par(gap.after = c(rep(5, length(unique(df[[1]]))-1), 15, 
                                   rep(5, length(unique(df[[2]]))-1), 15))
      }
      
      row = rep(brewer.pal(12, "Paired"), nrow(df))[1:length(unique(df[[1]]))]
      col = rep(rev(brewer.pal(12, "Paired")), nrow(df))[1:length(unique(df[[2]]))]
      
      # plot!
      chordDiagram(df, annotationTrack="grid", preAllocateTracks = list(track.height=0.2), grid.col=c(row,col))
      title(sampleName, cex=0.8)
      circos.trackPlotRegion(track.index = 1, bg.border = NA,
                             panel.fun = function(x, y) {
                               sector.name = get.cell.meta.data("sector.index")
                               xlim = get.cell.meta.data("xlim")
                               ylim = get.cell.meta.data("ylim")
                               circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", adj = c(0, 1.5))
                             }
      )
      circos.clear()
      dev.off()
  }
}

plotCirclize("LambdaR2_L001", "/Users/u0001382/sandbox/LambdaR2_BGPV9_AGGCATCT-ATCACG_L001/abundance/")
