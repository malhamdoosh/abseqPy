#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: Plots a circos diagram of V-J associations                         #
#################################################################################

library(circlize)

circulize <- function(sample_name, dir, path_to_clone) {
  
  file <- paste(dir, path_to_clone, sample_name, "_vjassoc.csv", sep="")
  
  df <- read.csv(file)
  
  circos.par(gap.after = c(rep(5, length(unique(df[[1]]))-1), 15, 
                           rep(5, length(unique(df[[2]]))-1), 15))
  
  row = rep(brewer.pal(12, "Paired"), nrow(df))[1:length(unique(df[[1]]))]
  col = rep(brewer.pal(12, "Paired"), nrow(df))[1:length(unique(df[[2]]))]
  chordDiagram(df, annotationTrack="grid", preAllocateTracks = list(track.height=0.2), grid.col=c(row,col))
  title(sample_name, cex=0.8)
  
  circos.trackPlotRegion(track.index = 1, bg.border = NA,
                         panel.fun = function(x, y) {
                           sector.name = get.cell.meta.data("sector.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", adj = c(0, 1.5))
                         }
  )
  
  circos.clear()
}

circos_plotter <- function(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone) {
    layout(matrix(1:3, 1, 3))
    circulize(sample_name1, dir1, path_to_clone)
    circulize(sample_name2, dir2, path_to_clone)
    circulize(sample_name3, dir3, path_to_clone)
}
