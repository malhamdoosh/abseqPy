#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: plots cdr3 spectratype of 3 samples in 'dodging' fashion           #
#################################################################################

library(splines)
library(ggplot2)


spectratype_plot <- function(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone) {
    file1 <- paste(dir1, path_to_clone, sample_name1, "_cdr3_spec.csv", sep="")
    file2 <- paste(dir2, path_to_clone, sample_name2, "_cdr3_spec.csv", sep="")
    file3 <- paste(dir3, path_to_clone, sample_name3, "_cdr3_spec.csv", sep="")

    df1 <- read.csv(file1)
    df2 <- read.csv(file2)
    df3 <- read.csv(file3)

    df1$round <- rep(sample_name1, nrow(df1))
    df2$round <- rep(sample_name2, nrow(df2))
    df3$round <- rep(sample_name3, nrow(df3))
    total1 <- sum(df1$count)
    total2 <- sum(df2$count)
    total3 <- sum(df3$count)
    df1$percent <- df1$count/total1
    df2$percent <- df2$count/total2
    df3$percent <- df3$count/total3

    df.inter_union <- rbind(df1, df2)
    df.union <- rbind(df.inter_union, df3)

    g <- ggplot(df.union, aes(length, percent)) +
      geom_bar(stat="identity", aes(fill=round), width=0.5, position="dodge") +
      #geom_smooth(aes(colour=round), se=F, method="glm", formula=y~ns(x, 3), lwd=0.7)+
      #geom_text_repel(aes(label = count), size = 3) +
      #geom_text(aes(label = count), vjust=-1, size = 3) +
      labs(title="CDR-3 amino acid spectratype",
           subtitle="Distribution of CDR-3 amino acid lengths",
           x="Length(AA)",
           y="Distribution")
    plot(g)
}
