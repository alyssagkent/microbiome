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
setwd("/workdir/users/agk85/CDC2/bins/diversity_figures")
infile = paste(genetype, "_", minreads,"_connection_counts.txt",sep="")
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
df$Connectivity = df$Connections/(df$Taxa*df$Genes)
df$Patient = sapply(strsplit(as.character(df$Sample),"-"), function(x) x[1])
df2 <- subset(df, Level =="s__")

#diversity from metaphlan

metaphlan<- read.table("/workdir/users/agk85/CDC2/metaphlan/cdc/mgm/CDC_mgm_metaphlan.txt",stringsAsFactors = FALSE,header=F,sep="\t")
header = metaphlan[1,]
colnames(metaphlan)<-as.character(unlist(header))
metaphlan <- metaphlan[-1,]
taxnames = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
levs = c(2,3,4,5,6,7)

#this gets a single name from the ID
l = 7
taxalev=taxnames[l]
levels <-as.vector(sapply(strsplit(as.character(metaphlan$ID),"\\|"), function(x) length(x)))
metaphlan$levels <- levels
metaphlan_sub <- subset(metaphlan, levels == l)

taxa <-sapply(strsplit(as.character(metaphlan_sub$ID),"__"), `[`, length(strsplit(as.character(metaphlan_sub$ID),"__")[[l+1]]))
rownames(metaphlan_sub)<-taxa

metaphlan_sub$levels <-NULL
metaphlan_sub$ID <-NULL
Patient<-sapply(strsplit(colnames(metaphlan_sub),"\\."), function(x) x[1])
patients <- factor(Patient)
pat <- data.frame(patients)

dft <- t(metaphlan_sub)
dflev <- matrix(as.numeric(unlist(dft)),nrow=nrow(dft))
colnames(dflev) <- colnames(dft)
rownames(dflev) <- rownames(dft)
#dflev$pat <- factor(Patient)
dflev.diversity <- diversity(dflev, index="shannon")
dflev.div.mat <- data.frame(names(dflev.diversity), dflev.diversity)
colnames(dflev.div.mat)<-c("Sample","Shannon_diversity")
mergedf <- merge(dflev.div.mat, df2, by=c("Sample"))

#get the number of metaphlan taxa
#make a boolean 
metaphlan_bool <- data.matrix(metaphlan_sub)
metaphlan_bool[metaphlan_bool>0]<-1
metaphlan_ntaxa <- data.frame(names(colSums(metaphlan_bool)),colSums(metaphlan_bool))
colnames(metaphlan_ntaxa)<-c("Sample","NTaxa")
mergedf2 <- merge(mergedf, metaphlan_ntaxa,by=c("Sample"))

#enteros
l = 5
taxalev=taxnames[l]
levels <-as.vector(sapply(strsplit(as.character(metaphlan$ID),"\\|"), function(x) length(x)))
metaphlan$levels <- levels
enteros <- subset(metaphlan, ID == "k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae")
enteros$levels <-NULL
enteros$ID <-NULL

metaphlan_enteros<-data.frame(names(enteros),t(enteros))
colnames(metaphlan_enteros)<-c("Sample","Enterobacteriaceae")
mergedf3 <- merge(mergedf2, metaphlan_enteros, by=c("Sample"))

tps <- read.csv("/workdir/users/agk85/CDC2/MetaDesignTPs.txt", header=F, sep="\t")
colnames(tps)<- c("Sample","TP")

readcounts = read.table("/workdir/users/agk85/CDC2/read_distributions/ReadCounts.txt",sep="\t",header=T)


############################################33
#nine panel figures
#xaxis timepoints 
#1 shannon
#2 num taxa
#3 num species metaphlan
#4 connectivity

df3tp <- merge(mergedf3, tps, by="Sample")
df3tprc <- merge(df3tp, readcounts, by="Sample")

#import ARG RPKM
infile2=read.csv("/workdir/users/agk85/CDC2/args/mapping/bwa_alignments_99_99/arg_v_samp_99_99.txt",header =T,row.names=1)
total_rpkm <- rowSums(infile2)
tr <- total_rpkm[order(names(total_rpkm))]
tot_rpkm <- data.frame(names(tr),tr)
colnames(tot_rpkm)<-c("Sample","Total_ARG_RPKM")

df3tprcta <- merge(df3tprc,tot_rpkm,by=c("Sample"))


df3tprcta$Level <- NULL
df.melt <- melt(df3tprcta, id.vars=c("Sample","Patient","TP"))

outhandle = paste("Allmetrics_taxa_",genetype,"_",minreads,".pdf",sep="")
pdf(outhandle,height=20,width=15,useDingbats=FALSE)
ggplot(df.melt,aes(x=as.factor(TP),y=as.numeric(as.character(value)),group=Patient))+
geom_point()+
geom_line(aes(col=Patient))+
facet_grid(variable~Patient,scales="free")+
scale_color_manual(values=kewl)+
xlab("Timepoints")+
ylab("")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



#regression
library (lmerTest) # Mixed model package by Douglas Bates, comes w/ pvalues! 
library (texreg) #Helps us make tables of the mixed models
library (afex) # Easy ANOVA package to compare model fits
library (plyr) # Data manipulator package
library (ggplot2) # GGplot package for visualizing data

lm_eqn <- function(df){
    m <- lm(Connectivity ~ NTaxa, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

df<- subset(df3tprcta, select=c("Ntaxa","Connectivity"))
pdf("NTaxa_v_connectivity.pdf",height=5,width=8,useDingbats=F)
ggplot(df3tprcta, aes(x=NTaxa, y=Connectivity))+
geom_point(aes(color=Patient))+
scale_color_manual(values=kewl)+
xlab("Metaphlan Species")+
ylab("Connectivity")+
#geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
#geom_text(x=60,y=0.07, label=lm_eqn(df), parse = TRUE)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()




library(lme4)
library(nlme)
library (lmerTest)
#merge the sick healthy
sh <- cbind(sort(unique(df3tprcta$Patient)),c("Sick","Sick","Sick","Sick","Sick","Sick","Sick","Healthy","Healthy")
colnames(sh)<-c("Patient","Status")

df3tprcta <- merge(df3tprcta,sh,by=c("Patient"),all.x=T)
df3tprcta$TP<-as.factor(df3tprcta$TP)
model <- lmer(Connectivity ~ NTaxa+Status + (1|Patient/TP), data=df3tprcta)

m1 <- lmer(Connectivity ~ NTaxa+Status + (1|Patient) +(1|TP), data=df3tprcta)
library(lmerTest)
summary(m1)
anova(m1)


# pvalues for this are weird so...be careful


df3tprcta$fit <- predict(m1)


pdf("Boxplot_connectivity.pdf",height=5, width=8,useDingbats=F)
ggplot(df3tprcta,aes(x=Patient, y=Connectivity))+
geom_boxplot(outlier.shape = NA)+
geom_point(aes(color=Patient,position=position_jitter(0.2))+
scale_color_manual(values=kewl)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf("NTaxa_v_connectivity_fit.pdf",height=5,width=8,useDingbats=F)
ggplot(df3tprcta, aes())+
geom_point(aes(x=NTaxa,y=Connectivity,color=Patient))+
scale_color_manual(values=kewl)+
geom_line(aes(x=NTaxa,y=fit,col=Patient,group=Patient))+
xlab("Metaphlan Species")+
ylab("Connectivity")+
#geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
#geom_text(x=60,y=0.07, label=lm_eqn(df), parse = TRUE)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()











