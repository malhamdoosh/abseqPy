#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: calculates the cardinality of intersecting clones across           #
#            3 or 2 samples                                                     #
#################################################################################

library(VennDiagram)
library(gridExtra)

venn2samples <- function(sample_name1, sample_name2, dir1, dir2, path_to_clone, top) {
    file1 <- paste(dir1, path_to_clone, sample_name1, "_cdr3_clonotypes_100_over.csv.gz", sep="")
    file2 <- paste(dir2, path_to_clone, sample_name2, "_cdr3_clonotypes_100_over.csv.gz", sep="")
    df1 <- read.csv(file1, stringsAsFactors=FALSE)
    df2 <- read.csv(file2, stringsAsFactors=FALSE)

    # some cleanup
    keep_col_name <- c("Clonotype", "Count")
    df1 <- df1[keep_col_name]
    df2 <- df2[keep_col_name]

    if (!missing(top)) {
        df1 <- head(df1, top)
        df2 <- head(df2, top)
    }

    # gather all clonotypes from repertoires
    df.union <- merge(df1, df2, by="Clonotype", all.y=TRUE, all.x=TRUE)
    colnames(df.union) <- c("Clonotype", sample_name1, sample_name2)
    # replace NaN with 0
    df.union[is.na(df.union)] <- 0

    countsample1 <- sum(df.union[sample_name1] > 0)
    countsample2 <- sum(df.union[sample_name2] > 0)


    countinter <- sum(df.union[sample_name1] > 0 & df.union[sample_name2] > 0)

    p <- draw.pairwise.venn(countsample1, countsample2, countinter,
      category = c(sample_name1, sample_name2),
      lty="blank",
      col="transparent",
      cat.fontface = "bold",
      cex=1.7, fill=2:3 , alpha=0.5,
      ind=FALSE,
      filename="something.pdf")
    return (p)
}

venn3samples <- function(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone, top) {

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

    if (!missing(top)) {
        df1 <- head(df1, top)
        df2 <- head(df2, top)
        df3 <- head(df3, top)
    }

    # gather all clonotypes from repertoires
    df.halfunion <- merge(df1, df2, by="Clonotype", all.y=TRUE, all.x=TRUE)
    colnames(df.halfunion) <- c("Clonotype", sample_name1, sample_name2)
    df.union <- merge(df.halfunion, df3, by="Clonotype", all.y=TRUE, all.x=TRUE)
    colnames(df.union) <- c("Clonotype", sample_name1, sample_name2, sample_name3)
    # replace NaN with 0
    df.union[is.na(df.union)] <- 0

    countsample1 <- sum(df.union[sample_name1] > 0)
    countsample2 <- sum(df.union[sample_name2] > 0)
    countsample3 <- sum(df.union[sample_name3] > 0)

    count1inter2 <- sum(df.union[sample_name1] > 0 & df.union[sample_name2] > 0)
    count1inter3 <- sum(df.union[sample_name1] > 0 & df.union[sample_name3] > 0)
    count2inter3 <- sum(df.union[sample_name2] > 0 & df.union[sample_name3] > 0)

    countinter <- sum(df.union[sample_name1] > 0 & df.union[sample_name2] > 0 & df.union[sample_name3] > 0)

    p <- draw.triple.venn(area1=(countsample1),
     area2=(countsample2),
     area3=(countsample3),
      n12=(count1inter2),
      n13=(count1inter3),
      n23=(count2inter3),
      n123=countinter,
      category = c(sample_name1, sample_name2, sample_name3),
      lty="blank",
      ind=FALSE,
      col="transparent", cat.fontface = "bold", cex=1.7, cat.cex=1.7, fill=2:4 , alpha=0.5, filename="clonotypeIntersection.pdf")
    return (p)
}

### USAGE ###
# p.all <- venn2diagram(file1, file2)
# grid.arrange(gTree(children=p.top), gTree(children=p.all), ncol=2, top=textGrob("Top 101 and all (CDR3)clones intersection", gp=gpar(fontsize=20,font=8)))
