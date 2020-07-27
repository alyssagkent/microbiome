#plot histograms for arg org and org arg relationships
#updated to only include things that have taxonomy down to species to be counted
library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(scales)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC2/arg_v_org/metagenomes/comparisons/nocorecorenomobilemobile/")

links = c(0,1,2)
for (link in links){

infile = paste('arg_orgcounts_99_0_species', link, '.txt', sep="")
indata = read.csv(infile, header=T, sep="\t")

#remove sample and Taxa
ind <- subset(indata, select=c(Cluster, ARGname, Study, species_count))
aggdata = aggregate(species_count ~ ., FUN=mean, na.rm=TRUE, data = ind)

pdf(paste("Healthy_CDC_averaged_boxplots_", link, "link.pdf", sep = ""))
ggplot(aggdata,aes(x = Study, y=species_count)) +
	geom_boxplot()+
	geom_point(position = position_jitter())
dev.off()



aggcast <- cast(Cluster+ARGname ~ Study, data=aggdata)
t = t.test(x=aggcast$CDC, y=aggcast$Healthy, paired=T)
print(t)

}

df <- data.frame(group, session, value, index, U = interaction(session, group))
p <- ggplot(df, aes(x = U, y = value, fill = session)) + 
  scale_x_discrete(labels = rep(unique(group), each = 2))
p <- p + geom_line(aes(group = index), alpha = 0.6, colour = "black", data = df) 
