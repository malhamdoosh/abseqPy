library(ggplot2)
source('~/Documents/repo/abseq/templates/summarySE.R')
theme_set(theme_bw())

sample_name1 <- "PCR1_L001"
sample_name2 <- "PCR2_L001"
sample_name3 <- "PCR3_L001"
dir1 <- "PCR1_BH5C6_CAACG-CGTGAT_L001"
dir2<- "PCR2_BH5C6_CTAGCT-CTAGTACG_L001"
dir3 <-"PCR3_BH5C6_AGGCAGAA-TACAGC_L001"
path_to_div <- "/diversity/"
f1 <- paste(dir1, path_to_div, sample_name1, "_cdr_v_recapture.csv.gz", sep="")
f2 <- paste(dir2, path_to_div, sample_name2, "_cdr_v_recapture.csv.gz", sep="")
f3 <- paste(dir3, path_to_div, sample_name3, "_cdr_v_recapture.csv.gz", sep="")

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
df.1 <- df.1[df.1$region %in% c("CDR3", "V"), ]
df.2 <- df.2[df.2$region %in% c("CDR3", "V"), ]
df.3 <- df.3[df.3$region %in% c("CDR3", "V"), ]

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
  labs(title = "Percent Recapture of CDRs and V Domains",
       subtitle="Mean number of recaptured sequences with 95% confidence interval",
       x = "Sample size", y = "Percent Recapture")
   
plot(p)