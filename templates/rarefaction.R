library(ggplot2)
theme_set(theme_bw())

# from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

sample_name1 <- "PCR1_L001"
sample_name2 <- "PCR2_L001"
sample_name3 <- "PCR3_L001"
dir1 <- "PCR1_BH5C6_CAACG-CGTGAT_L001"
dir2<- "PCR2_BH5C6_CTAGCT-CTAGTACG_L001"
dir3 <-"PCR3_BH5C6_AGGCAGAA-TACAGC_L001"
path_to_div <- "/diversity/"
f1 <- paste(dir1, path_to_div, sample_name1, "_cdr_v_rarefaction.csv.gz", sep="")
f2 <- paste(dir2, path_to_div, sample_name2, "_cdr_v_rarefaction.csv.gz", sep="")
f3 <- paste(dir3, path_to_div, sample_name3, "_cdr_v_rarefaction.csv.gz", sep="")

# first, obtain the x-axis labels and ticks (metadata from csv's first row) - infact, we wan't the maximum (temporary solution for now)
# things are going to get ugly when the sample sizes are vastly different (you can barely see the tiny samples)
fp <- file(f1, "r")
xticks <- eval(parse(text=readLines(fp, n=1)))
close(fp)

fp <- file(f2, "r")
candidate <- eval(parse(text=readLines(fp, n=1)))
close(fp)
if (tail(candidate, n=1) > tail(xticks, n=1)) {
  xticks <- candidate
}

fp <- file(f3, "r")
candidate <- eval(parse(text=readLines(fp, n=1)))
close(fp)
if (tail(candidate, n=1) > tail(xticks, n=1)) {
  xticks <- candidate
}
# done

# read files into dataframes
df.1 <- read.csv(f1, skip=1)
df.2 <- read.csv(f2, skip=1)
df.3 <- read.csv(f3, skip=1)


df.1$sample <- rep(sample_name1, nrow(df.1))
df.2$sample <- rep(sample_name2, nrow(df.2))
df.3$sample <- rep(sample_name3, nrow(df.3))

# uncomment to get CDR3 and V region only (by default, all CDR regions)
# df.1 <- df.1[df.1$region %in% c("CDR3", "V"), ]
# df.2 <- df.2[df.2$region %in% c("CDR3", "V"), ]
# df.3 <- df.3[df.3$region %in% c("CDR3", "V"), ]

# get mean, sd, se, and ci
df.1.m <- summarySE(df.1, measurevar='y', groupvars = c("x", "region", "sample"))
df.2.m <- summarySE(df.2, measurevar='y', groupvars = c("x", "region", "sample"))
df.3.m <- summarySE(df.3, measurevar='y', groupvars = c("x", "region", "sample"))


df.al <- rbind(df.1.m, df.2.m)
df.all <- rbind(df.al, df.3.m)
 
p <- ggplot(df.all, aes(x=x, y=y)) +
  geom_line(aes(linetype=region, color=sample)) +
  scale_x_continuous(breaks=xticks) +
  geom_ribbon(aes(ymin=y-ci, ymax=y+ci, fill=region), alpha=0.2) +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  labs(title = "Rarefaction of CDRs and V Domains",
       subtitle="Plot of mean number of deduplicated sequences with 95% confidence interval",
       x = "Sample size", y = "Number of Deduplicated Sequences")
   
plot(p)