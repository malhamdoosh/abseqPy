#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: Plots a stacked bar graph showing clonotype distributions side-by- #
#            side with other samples.                                           #
#################################################################################

library(ggplot2)
library(RColorBrewer)


litmus_paper <- function(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone) {
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

    # take top N rows
    df1 <- head(df1, 10)
    df2 <- head(df2, 10)
    df3 <- head(df3, 10)

    df1$round <- rep(sample_name1, nrow(df1))
    df2$round <- rep(sample_name2, nrow(df2))
    df3$round <- rep(sample_name3, nrow(df3))

    # get total number of clonotypes
    pcr1_count <- sum(df1$Count)
    pcr2_count <- sum(df2$Count)
    pcr3_count <- sum(df3$Count)

    # normalize values in table
    df1$Count <- df1$Count/pcr1_count
    df2$Count <- df2$Count/pcr2_count
    df3$Count <- df3$Count/pcr3_count

    df.interunion <- rbind(df1, df2)
    df.union <- rbind(df.interunion, df3)
    colourCount <- length(unique(df.union$Clonotype))
    #
    # taken from https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes/41230685 and modified
    c30 <- c("dodgerblue2","#E31A1C", # red
                "green4",
                "#6A3D9A", # purple
                "#FF7F00", # orange
                "black","gold1",
                "skyblue2","#FB9A99", # lt pink
                "olivedrab1",
                "#CAB2D6", # lt purple
                "#FDBF6F", # lt orange
                "gray70", "khaki2",
                "maroon","orchid1","deeppink1","blue1","steelblue4",
                "darkturquoise","green1","yellow4","yellow3",
                "aquamarine", "darkorange4", "mediumpurple1", "dimgrey", "darkseagreen1", "lightyellow", "coral2")
      
    g <- ggplot(df.union, aes(x=round, y=Count)) +
      geom_bar(stat='identity', aes(fill=Clonotype)) +
      theme(legend.position="bottom", legend.box = "horizontal", legend.title=element_blank(), legend.text=element_text(size=9)) +
      labs(title="Top 10 clonotype across each sample",
           subtitle="Colour coded clonotypes, distribution of each clonotype is relative to top10, not overall.",
           x="round",
           y="Distribution")  +
    scale_fill_manual(values=c30)

    plot(g)
}
