require(reshape2)
require(ggplot2)

file = "/Users/harry/sandbox/small/test/results/PCR3_BH5C6_AGGCAGAA-TACAGC_L001/diversity/PCR3_L001_clonotype_diversity_region_analysis.csv.gz"
df <- read.csv(file, stringsAsFactors=FALSE)

headers = c("fr1", "cdr1", "fr2", "cdr2", "fr3", "fr4")

# sort the df with decreasing counts of CDR3 occurance
df <- df[with(df, order(-count)), ]

# grab top N
df <- head(df, 15)

# add new column to sum the "unique" regions
df$sumcounts = rowSums(df[, headers])

# add new row as reference of "unique regions"
de <- data.frame("", Inf,1,1,1,1,1,1,6)
names(de) <- names(df)
df <- rbind(df, de)


# reshape to multi-row
df.mel <- melt(df, measure.vars = headers)
df.mel[, "value"] = df.mel[, "value"] / df.mel[, "sumcounts"]

# plot!
g <- ggplot(df.mel, aes(x = cdr3, y = value, fill=variable, label = sprintf("%0.2f%%", round(value*100, digits=2)))) +
  geom_bar(stat='identity') +
  theme(text=element_text(size=10), axis.text.x = element_text(angle=65, hjust=1)) +
  #geom_text(position = position_stack(vjust = 0.5)) +
  #stat_summary(fun.y = sum, aes(label = sumcounts, group=cdr3), geom='text', vjust=-.2) +
  labs(title="Varying levels of FRs and CDRs of top 15 CDR3 clonotype",
       subtitle="Counts of unique regions for a given CDR3",
       x="CDR3",
       y="Proportion") +
  guides(fill=guide_legend(title="Region"))

plot(g)