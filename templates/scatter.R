#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: plots a scatter plot of clonotypes in both samples.                #
#            (x, y) represents the cardinality of the particular clonotype in   #
#            repertoire (X) and repertoire (Y)                                  #
#################################################################################

library(gridExtra)
library(ggplot2)

scatter <- function(df1, df2, name1, name2, cloneClass) {
  # gather all clonotypes from both repertoire
  all_clones <- unique(c(df1$Clonotype, df2$Clonotype))
  df.union <- merge(df1, df2, by="Clonotype", all.y=TRUE, all.x=TRUE)
  # replace NaN with 0
  df.union[is.na(df.union)] <- 0
  
  # BEGIN PLOT
  options(scipen=999)
  library(ggplot2)
  theme_set(theme_bw())
  
  gg <- ggplot(df.union, aes(x=Count.x, y=Count.y)) +
    geom_point() +
    labs(subtitle=paste(name2, " vs ", name1, " plot based on clonotype counts"),
         y=name2,
         x=name1,
         title=paste("Scatter plot of ", cloneClass, " clonotype counts"))
  return(gg)
}

scatter_clones <- function(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone, cloneClass) {
    file1 <- paste(dir1, path_to_clone, sample_name1, "_cdr3_clonotypes_100_over.csv.gz", sep="")
    file2 <- paste(dir2, path_to_clone, sample_name2, "_cdr3_clonotypes_100_over.csv.gz", sep="")
    file3 <- paste(dir3, path_to_clone, sample_name3, "_cdr3_clonotypes_100_over.csv.gz", sep="")
    df1 <- read.csv(file1, stringsAsFactors=FALSE)
    df2 <- read.csv(file2, stringsAsFactors=FALSE)
    df3 <- read.csv(file3, stringsAsFactors=FALSE)

    # some cleanup
    keep_col_name <- c("Clonotype", "Count")
    df1 <- df1[keep_col_name]
    df2 <- df2[keep_col_name]
    df3 <- df3[keep_col_name]
    p1 <- scatter(df1, df2, sample_name1, sample_name2, cloneClass)
    #p2 <- scatter(df1, df3, sample_name1, sample_name3, cloneClass)
    p3 <- scatter(df2, df3, sample_name2, sample_name3, cloneClass)
    grid.arrange(p1, p3, ncol=2)
    #grid.arrange(p1, p2, p3, ncol=3)
}

