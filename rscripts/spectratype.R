library(ggplot2)

plotSpectratype <- function(dataframes, sampleNames, region = "CDR3") {
  # Plots an amino acid spectratype diagram for a given region.
  # Args:
  #     dataframes: A list() type. List of dataframes.
  #     sampleNames: A vector type. Vectors of strings each representing the samplename
  #     region: String that will be displayed in the plot title. This specifies which region this spectratype belongs to
  nsample <- length(dataframes)
  stopifnot(nsample == length(sampleNames))
  
  # pre-processing
  for (i in 1:nsample) {
    df <- dataframes[[i]]
    df$sample<- rep(sampleNames[i], nrow(df))
    df$percent <- df$count / sum(df$count)
    dataframes[[i]] <- df
  }
  
  # merge
  if (nsample > 1) {
    df.union <- rbind(dataframes[[1]], dataframes[[2]])
    if (nsample > 2) {
      for (i in 3:nsample) {
        df.union <- rbind(df.union, dataframes[[i]])
      }
    }
  } else {
    df.union <- dataframes[[1]]
  }
  
  g <- ggplot(df.union, aes(length, percent)) +
    geom_bar(stat = "identity", aes(fill = sample), width = 0.5, position = "dodge") +
    #geom_smooth(aes(colour=round), se=F, method="glm", formula=y~ns(x, 3), lwd=0.7)+
    #geom_text_repel(aes(label = count), size = 3) +
    #geom_text(aes(label = count), vjust=-1, size = 3) +
    labs(title = paste(region, "amino acid spectratype"),
         subtitle = paste("Distribution of", region, "amino acid lengths"),
         x = "Length(AA)",
         y = "Distribution")
  return (g)
}

# fs <- list.files(pattern = "PCR[123].*_cdr3_spectratype_no.*.csv", recursive = TRUE, full.names = TRUE)
# p <- plotSpectratype(lapply(fs, read.csv, stringsAsFactors = FALSE), c("PCR1", "PCR2", "PCR3"), "CDR3")
# plot(p)