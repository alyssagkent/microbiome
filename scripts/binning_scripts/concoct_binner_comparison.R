#plot histograms for arg org and org arg relationships
#updated to only include things that have taxonomy down to species to be counted
#cat *-1/checkm_lineage/*.named > concoct_t1_checkm_stats.txt
#cat *-1/checkm_lineage/*.named > concoct_checkm_stats.txt
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
setwd("/workdir/users/agk85/CDC2/concoct")

header = c("Bin","Taxa","Completeness","Contamination","Hetero","Size")

infile = read.csv("concoct_t1_checkm_stats.txt", sep="\t",header=F)
colnames(infile)<-header
infile$Binner = rep("Concoct_t1",length(infile$Bin))

infile1 = read.csv("concoct_checkm_stats.txt", sep="\t", header =F)
colnames(infile1)<-header
#infile1$Sample = as.vector(sapply(strsplit(as.character(infile$Bin), "_"), function(x) x[[1]]))
infile1$Binner = rep("Concoct", length(infile1$Bin))

infiles = rbind(infile, infile1)
dim(infiles)

df <- subset(infiles, Size>100000)
dim(df)

pdf("Genome_bin_sizes4.pdf")
ggplot(df, aes(Binner, Size))+
geom_boxplot(outlier.colour=NA)+
geom_point(position=position_jitter(width=.2))+
scale_y_log10(labels = comma)
dev.off()


pdf("Genome_bin_contamination4.pdf")
ggplot(df, aes(Binner, Contamination))+
geom_boxplot(outlier.colour=NA)+
geom_point(position=position_jitter(width=.2))+
scale_y_continuous(labels = comma)
dev.off()

pdf("Genome_bin_contamination4_0_200.pdf")
ggplot(df, aes(Binner, Contamination))+
geom_boxplot(outlier.colour=NA)+
geom_point(position=position_jitter(width=.2))+
scale_y_continuous(labels = comma,limits=c(0,200))
dev.off()


pdf("Genome_bin_completeness4.pdf")
ggplot(df, aes(Binner, Completeness))+
geom_violin()+
geom_point(position=position_jitter(width=.1))+
scale_y_continuous(labels = comma)
dev.off()

library(psych)
a = describeBy(df, "Binner")

hq <- subset(df, Contamination<5 & Completeness>90)
mq <- subset(df, Contamination<10 & Completeness>=50)
lq <- subset(df, Contamination<10 & Completeness<50)
bad<-subset(df, Contamination>=10)

b = describeBy(hq, "Binner")
c = describeBy(mq, "Binner")
d = describeBy(lq, "Binner")
e = describeBy(bad, "Binner")


hqstats = rbind(b$Concoct, b$Concoct_t1)
mqstats = rbind(c$Concoct, c$Concoct_t1)
lqstats = rbind(d$Concoct, d$Concoct_t1)

write.csv(hqstats, "HQ_bins_stats.txt")
write.csv(mqstats, "MQ_bins_stats.txt")
write.csv(lqstats, "LQ_bins_stats.txt")
