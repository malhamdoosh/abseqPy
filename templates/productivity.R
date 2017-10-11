#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: plots a facet_grid plot to show n-samples with productivity        #
#            comparision (and reason - stop/frame/both)                         #
#################################################################################

library(ggplot2)
library(ggrepel)

productivty_plot <- function(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone) {

    file1 <- paste(dir1, path_to_clone, sample_name1, "_productivity2.csv", sep="")
    file2 <- paste(dir2, path_to_clone, sample_name2, "_productivity2.csv", sep="")
    file3 <- paste(dir3, path_to_clone, sample_name3, "_productivity2.csv", sep="")

    df1 <- read.csv(file1)
    df2 <- read.csv(file2)
    df3 <- read.csv(file3)

    df1$round <- rep(sample_name1, nrow(df1))
    df2$round <- rep(sample_name2, nrow(df2))
    df3$round <- rep(sample_name3, nrow(df3))

    df.inter_union <- rbind(df1, df2)
    df.union <- rbind(df.inter_union, df3)
    drops <- c("X")
    df.union <- df.union[, !(names(df.union) %in% drops)]

    g <- ggplot(df.union, aes(round, Percentage)) +
      geom_bar(stat="identity", aes(fill=Reason), width=0.5) +
      facet_grid(~ Productivity)+
      labs(title="Productivity",
          subtitle="Percentage of unproductive reads due to stop codons and frameshifts",
          x="Round",
          y="Percentage") +
      scale_y_continuous(limits=c(0,100))
    plot(g)
}
