library(ggplot2)

theme_set(theme_bw())
sample_name1 <- "PCR1_L001"
sample_name2 <- "PCR2_L001"
sample_name3 <- "PCR3_L001"
dir1 <- "PCR1_BH5C6_CAACG-CGTGAT_L001"
dir2<- "PCR2_BH5C6_CTAGCT-CTAGTACG_L001"
dir3 <-"PCR3_BH5C6_AGGCAGAA-TACAGC_L001"
path_to_div <- "/diversity/"
f1 <- paste(dir1, path_to_div, sample_name1, "_cdr_v_duplication.csv", sep="")
f2 <- paste(dir2, path_to_div, sample_name2, "_cdr_v_duplication.csv", sep="")
f3 <- paste(dir3, path_to_div, sample_name3, "_cdr_v_duplication.csv", sep="")

# first, obtain the x-axis labels and ticks (metadata from csv's first 2 rows)
fp <- file(f1, "r")
xticks <- eval(parse(text=readLines(fp, n=1)))
xlabels <- eval(parse(text=readLines(fp, n=1)))
close(fp)

# read files into dataframes
df.1 <- read.csv(f1, skip=2)
df.2 <- read.csv(f2, skip=2)
df.3 <- read.csv(f3, skip=2)


df.1$sample <- rep(sample_name1, nrow(df.1))
df.2$sample <- rep(sample_name2, nrow(df.2))
df.3$sample <- rep(sample_name3, nrow(df.3))

# uncomment to get CDR3 and V region only (by default, all CDR regions)
df.1 <- df.1[df.1$region %in% c("CDR3", "V"), ]
df.2 <- df.2[df.2$region %in% c("CDR3", "V"), ]
df.3 <- df.3[df.3$region %in% c("CDR3", "V"), ]


df.al <- rbind(df.1, df.2)
df.all <- rbind(df.al, df.3)

p <- ggplot(df.all, aes(x=x, y=y)) +
  geom_line(aes(linetype=region, color=sample)) +
  scale_x_continuous(breaks=xticks, labels=xlabels) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  labs(title = "Duplication levels of CDRs and V Domains", x = "Duplication Level", y = "Proportion of Duplicated Sequences")
  
plot(p)