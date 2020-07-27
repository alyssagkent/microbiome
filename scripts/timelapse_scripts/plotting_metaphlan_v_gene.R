library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(scales)
set.seed(42)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.grid = element_blank())}
theme_set(theme_nogrid())

minreads = 2
genetype = 'arg'

minreads=0
args <- commandArgs(trailingOnly = TRUE)
minreads = args[1]
genetype = args[2]

###########
#import data
###########
setwd("/workdir/users/agk85/CDC2/bins/gene_abundances")
infile = paste("metaphlan_v_generpkm_",genetype, "_org_", minreads,".txt",sep="")
df <-read.table(infile,header=T, sep="\t")
################
#colors
################
arg.palette <-colorRampPalette(brewer.pal(12,"Set3"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 9
kewl<-col_vector[60:(n+60)]
################################
df$level<- factor(df$level, levels = c('k__','p__','c__','o__','f__','g__', 's__'))


outhandle = paste("metaphlan_v_generpkm_",genetype,"_",minreads,".pdf",sep="")
pdf(outhandle,height=20,width=15,useDingbats=FALSE)
ggplot(df,aes(x=taxon_abund,y=generpkm, group=patient))+
geom_point(aes(col=patient))+
facet_grid(level~.,scales="free")+
scale_color_manual(values=kewl)+
xlab("Taxon relative abundance (Metaphlan)")+
ylab("Gene abundance (RPKM)")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


outhandle = paste("metaphlan_v_generpkm_",genetype,"_",minreads,"_0_200.pdf",sep="")
pdf(outhandle,height=20,width=15,useDingbats=FALSE)
ggplot(df,aes(x=taxon_abund,y=generpkm, group=patient))+
geom_point(aes(col=patient))+
facet_grid(level~.,scales="free")+
scale_color_manual(values=kewl)+
xlab("Taxon relative abundance (Metaphlan)")+
ylab("Gene abundance (RPKM)")+
ylim(0,200)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

outhandle = paste("metaphlan_v_generpkm_",genetype,"_",minreads,"_0_10.pdf",sep="")
pdf(outhandle,height=20,width=15,useDingbats=FALSE)
ggplot(df,aes(x=taxon_abund,y=generpkm, group=patient))+
geom_point(aes(col=patient))+
facet_grid(level~.,scales="free")+
scale_color_manual(values=kewl)+
xlab("Taxon relative abundance (Metaphlan)")+
ylab("Gene abundance (RPKM)")+
ylim(0,10)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()





df_no0 <- subset(df, generpkm>0)
outhandle = paste("histogram_generpkm_",genetype,"_",minreads,"_0_20.pdf",sep="")
pdf(outhandle,height=20,width=15,useDingbats=FALSE)
ggplot(df_no0,aes(x=generpkm, fill=patient))+
geom_histogram(binwidth=0.5)+
facet_grid(level~.,scales="free")+
xlim(c(0,20))+
xlab("Gene abundance (RPKM)")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

outhandle = paste("histogram_taxon_abund_",genetype,"_",minreads,"_0_100.pdf",sep="")
pdf(outhandle,height=10,width=5,useDingbats=FALSE)
ggplot(df_no0,aes(x=taxon_abund, fill=patient))+
geom_histogram(binwidth=0.5)+
facet_grid(level~.,scales="free")+
xlim(c(0,100))+
xlab("Taxon abundance (Metaphlan)")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

