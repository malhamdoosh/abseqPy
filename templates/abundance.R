#################################################################################
#   Author : JiaHong FONG                                                       #
#   Date   : Tue 26 Sep 2017 15:27:05 AEST                                      #
#   Purpose: plots a dodging bar plot of 3 samples (family V-D-J) distribution  #
#################################################################################

library(ggplot2)
library(ggrepel)
library(gridExtra)
theme_set(theme_classic())

plot_dist <- function(df1, df2, df3, family, sample_name1, sample_name2, sample_name3, title="family") {
  CUTOFF <- 15
  # ### CLEANUP #### #
  # Don't need the total row (remove it)
  df1 <- df1[!(df1$Germline.group == 'TOTAL'), ]
  df2 <- df2[!(df2$Germline.group == 'TOTAL'), ]
  df3 <- df3[!(df3$Germline.group == 'TOTAL'), ]
  caps <- ""
  
  if (title != "family") {
    if (nrow(df1) > CUTOFF) {
      caps <- "Cutoff at top 15"
      df1 <- head(df1, CUTOFF)
    }
    if (nrow(df2) > CUTOFF) {
      caps <- "Cutoff at top 15"
      df2 <- head(df2, CUTOFF)
    }
    if (nrow(df3) > CUTOFF) {
      caps <- "Cutoff at top 15"
      df3 <- head(df3, CUTOFF)
    }
  }
  # ### CLEANUP DONE #### #
  
  # sort descending order (not necessary, but sort anyway).
  df1 <- df1[with(df1, order(-Count)), ]
  df2 <- df2[with(df2, order(-Count)), ]
  df3 <- df3[with(df3, order(-Count)), ]
  
  # add a column specifying where this data came from (which sample)
  df1$round <- rep(sample_name1, nrow(df1))
  df2$round <- rep(sample_name2, nrow(df2))
  df3$round <- rep(sample_name3, nrow(df3))
  
  df.intermediate_union <- rbind(df1, df2)
  df.union <- rbind(df.intermediate_union, df3)
  
  g <- ggplot(df.union, aes(Germline.group, Percentage....)) +
    geom_bar(stat="identity", aes(fill=round), width=0.5, position="dodge") +
    theme(text=element_text(size=10), axis.text.x = element_text(angle=74, vjust=0.4)) +
    #geom_text_repel(aes(label = Count), size = 3) +
    #geom_text(aes(label = Count), vjust=-1, size = 3) +
    labs(title=paste("Histogram of ", family, " ", title, " distribution", sep=""),
         subtitle=paste("Percentage of IGH", family, " ", title," distribution across all rounds", sep=""),
         x=paste(family,"-Germline ", title , sep=""),
         y="Percentage", caption=caps)
  return(g)
}

# ====================================================== #
# Family level stacked bar plot (Histogram)
# ====================================================== #
abundance_plot <- function(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone) {
  # V-FAMILY
  file1.v <- paste(dir1, path_to_clone, sample_name1, "_igv_dist_family_level.csv", sep="")
  file2.v <- paste(dir2, path_to_clone, sample_name2, "_igv_dist_family_level.csv", sep="")
  file3.v <- paste(dir3, path_to_clone, sample_name3, "_igv_dist_family_level.csv", sep="")
  df1.v <- read.csv(file1.v, stringsAsFactors=FALSE)
  df2.v <- read.csv(file2.v, stringsAsFactors=FALSE)
  df3.v <- read.csv(file3.v, stringsAsFactors=FALSE)
  
  # D-FAMILY
  file1.d <- paste(dir1, path_to_clone, sample_name1, "_igd_dist_family_level.csv", sep="")
  file2.d <- paste(dir2, path_to_clone, sample_name2, "_igd_dist_family_level.csv", sep="")
  file3.d <- paste(dir3, path_to_clone, sample_name3, "_igd_dist_family_level.csv", sep="")
  df1.d <- read.csv(file1.d, stringsAsFactors=FALSE)
  df2.d <- read.csv(file2.d, stringsAsFactors=FALSE)
  df3.d <- read.csv(file3.d, stringsAsFactors=FALSE)
  
  # J-FAMILY
  file1.j <- paste(dir1, path_to_clone, sample_name1, "_igj_dist_family_level.csv", sep="")
  file2.j <- paste(dir2, path_to_clone, sample_name2, "_igj_dist_family_level.csv", sep="")
  file3.j <- paste(dir3, path_to_clone, sample_name3, "_igj_dist_family_level.csv", sep="")
  df1.j <- read.csv(file1.j, stringsAsFactors=FALSE)
  df2.j <- read.csv(file2.j, stringsAsFactors=FALSE)
  df3.j <- read.csv(file3.j, stringsAsFactors=FALSE)
  
  p1 <- plot_dist(df1.v, df2.v, df3.v, "V", sample_name1, sample_name2, sample_name3)
  p2 <- plot_dist(df1.d, df2.d, df3.d, "D", sample_name1, sample_name2, sample_name3)
  p3 <- plot_dist(df1.j, df2.j, df3.j, "J", sample_name1, sample_name2, sample_name3)
  
  # ====================================================== #
  # Gene level stacked bar plot (Histogram)
  # ====================================================== #
  
  # V-GENE
  file1.v <- paste(dir1, path_to_clone, sample_name1, "_igv_dist_gene_level.csv", sep="")
  file2.v <- paste(dir2, path_to_clone, sample_name2, "_igv_dist_gene_level.csv", sep="")
  file3.v <- paste(dir3, path_to_clone, sample_name3, "_igv_dist_gene_level.csv", sep="")
  df1.v <- read.csv(file1.v, stringsAsFactors=FALSE)
  df2.v <- read.csv(file2.v, stringsAsFactors=FALSE)
  df3.v <- read.csv(file3.v, stringsAsFactors=FALSE)
  
  # D-GENE
  file1.d <- paste(dir1, path_to_clone, sample_name1, "_igd_dist_gene_level.csv", sep="")
  file2.d <- paste(dir2, path_to_clone, sample_name2, "_igd_dist_gene_level.csv", sep="")
  file3.d <- paste(dir3, path_to_clone, sample_name3, "_igd_dist_gene_level.csv", sep="")
  df1.d <- read.csv(file1.d, stringsAsFactors=FALSE)
  df2.d <- read.csv(file2.d, stringsAsFactors=FALSE)
  df3.d <- read.csv(file3.d, stringsAsFactors=FALSE)
  
  # J-GENE
  file1.j <- paste(dir1, path_to_clone, sample_name1, "_igj_dist_variant_level.csv", sep="")
  file2.j <- paste(dir2, path_to_clone, sample_name2, "_igj_dist_variant_level.csv", sep="")
  file3.j <- paste(dir3, path_to_clone, sample_name3, "_igj_dist_variant_level.csv", sep="")
  df1.j <- read.csv(file1.j, stringsAsFactors=FALSE)
  df2.j <- read.csv(file2.j, stringsAsFactors=FALSE)
  df3.j <- read.csv(file3.j, stringsAsFactors=FALSE)
  
  p1.p <- plot_dist(df1.v, df2.v, df3.v, "V", sample_name1, sample_name2, sample_name3, "gene")
  p2.p <- plot_dist(df1.d, df2.d, df3.d, "D", sample_name1, sample_name2, sample_name3, "gene")
  ### REMEMBER, abudance_plot(J ==> title="variant" NOT GENE)
  p3.p <- plot_dist(df1.j, df2.j, df3.j, "J", sample_name1, sample_name2, sample_name3, "variant")
  
  grid.arrange(p1, p2, p3, p1.p, p2.p, p3.p, ncol=3)
}


abundance_plot(sample_name1, sample_name2, sample_name3, dir1, dir2, dir3, path_to_clone)