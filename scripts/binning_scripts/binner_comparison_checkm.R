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
setwd("/workdir/users/agk85/CDC2/maxbin")

header = c("Bin","Taxa","Completeness","Contamination","Hetero","Size")

infile = read.csv("maxbin_checkm_stats.txt", sep="\t",header=F)
colnames(infile)<-header
#infile$Sample = as.vector(sapply(strsplit(as.character(infile$Bin), "[.]"),function(x) x[[1]]))
infile$Binner = rep("Maxbin",length(infile$Bin))

infile1 = read.csv("concoct_checkm_stats.txt", sep="\t", header =F)
colnames(infile1)<-header
#infile1$Sample = as.vector(sapply(strsplit(as.character(infile$Bin), "_"), function(x) x[[1]]))
infile1$Binner = rep("Concoct", length(infile1$Bin))

infile2 = read.csv("metabat_checkm_stats.txt", sep="\t",header=F)
colnames(infile2)<-header
#infile2$Sample = as.vector(sapply(strsplit(as.character(infile2$Bin), "_"),function(x) x[[1]]))
infile2$Binner = rep("Metabat",length(infile2$Bin))


infile3= read.csv("das_checkm_stats.txt", sep="\t", header=F)
colnames(infile3) <-header
#infile3$partial = as.vector(sapply(strsplit(as.character(infile3$Bin), "_"),function(x) x[[1]]))
#infile3$Sample = as.vector(sapply(strsplit(as.character(infile3$partial), "[.]"),function(x) x[[1]]))
#infile3$partial <- NULL
infile3$Binner = rep("DAS",length(infile3$Bin))

infile4= read.csv("das_withconcoct_checkm_stats.txt", sep="\t", header=F)
colnames(infile4) <-header
infile4$Binner = rep("DAS_with_Concoct",length(infile4$Bin))


infiles = rbind(infile, infile1, infile2, infile3,infile4)
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


hqstats = rbind(b$DAS_with_Concoct, b$DAS, b$Maxbin, b$Metabat, b$Concoct)
mqstats = rbind(c$DAS_with_Concoct, c$DAS, c$Maxbin, c$Metabat, c$Concoct)
lqstats = rbind(d$DAS_with_Concoct, d$DAS, d$Maxbin, d$Metabat, d$Concoct)

write.csv(hqstats, "HQ_bins_stats.txt")
write.csv(mqstats, "MQ_bins_stats.txt")
write.csv(lqstats, "LQ_bins_stats.txt")
