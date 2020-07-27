library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(pastecs)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.grid = element_blank())}
theme_set(theme_nogrid())
makeTransparent<-function(someColor, alpha=100)
{  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})}

setwd("/workdir/users/agk85/CDC2/read_distributions/")
indata = read.csv("All_length_hichits_mobility_abundance.txt",header = T,sep="\t")
#colnames(indata)<- c("Scfid","Sample","Length","Hic_hits","Mobile")
indata$Mobility <- as.character(indata$Mobility)
#subdata <- subset(indata, Sample=="B314-1")
colnames(indata)<-c("Scfid","Sample", "Length","Hic_hits", "Mobility", "RPKM", "Coverage", "Cis_hic_hits")

#bins
bins = read.table("/workdir/users/agk85/CDC2/das/all_DASTool_scaffolds2bin.txt",header=F, sep="\t")
colnames(bins) <- c("Scfid","Binid")
longname <- sapply(strsplit(as.character(bins$Binid),'_'), function(x) x[[1]])
bins$Sample <- sapply(strsplit(as.character(longname),'[.]'), function(x) x[[1]])

#kraken
bintable = read.table("/workdir/users/agk85/CDC2/das/all_bintables.txt",header=F, sep= "\t")
colnames(bintable) <- c("Binid","Patient","Sample","Quality","Taxa")

df1 <- merge(indata, bins, by="Scfid",all.y=T)
df2 <- merge(df1, bintable, by="Binid",all.y=T)

#within sample Bin representation
df3 <- subset(df2, select = c("Sample","Patient","Binid","Hic_hits","Cis_hic_hits"))
df3.agg <- aggregate(cbind(Hic_hits, Cis_hic_hits) ~ Binid+Sample+Patient, data=df3, sum)
df3.agg$one<-(df3.agg$Hic_hits>0)
df3.agg$two<-(df3.agg$Hic_hits>1)
df3.agg$five<-(df3.agg$Hic_hits>4)
df3.agg$zero<-(df3.agg$Hic_hits>-1)
df3.agg2 <- aggregate(cbind(zero,one,two,five) ~ Sample+Patient, data=df3.agg,sum)
df3.agg2$oneperc <- df3.agg2$one/df3.agg2$zero
df3.agg2$twoperc <- df3.agg2$two/df3.agg2$zero
df3.agg2$fiveperc <- df3.agg2$five/df3.agg2$zero
write.table(df3.agg2,row.names=FALSE,quote = F, "Bin_hic_coverage.txt",sep="\t")
df3.agg2$zero <-NULL
df3.agg2$one <- NULL
df3.agg2$two <- NULL
df3.agg2$five <- NULL
colnames(df3.agg2)<-c("Sample","Patient","One+","Two+","Five+")
write.table(stat.desc(df3.agg2),row.names=FALSE,quote = F,"Bin_hic_coverage_statistics.txt",sep="\t")
df3.agg2.melt <- melt(df3.agg2)
pdf("Bin_hic_coverage.pdf")
ggplot(df3.agg2.melt, aes(x=Patient, y=value))+
geom_boxplot(outlier.shape = NA)+
geom_point()+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
facet_grid(.~variable)
dev.off()

#within sample taxa representation
df3 <- subset(df2, select = c("Sample","Patient","Taxa","Hic_hits","Cis_hic_hits"))
df3.agg <- aggregate(cbind(Hic_hits, Cis_hic_hits) ~ Taxa+Sample+Patient, data=df3, sum)
df3.agg$one<-(df3.agg$Hic_hits>0)
df3.agg$two<-(df3.agg$Hic_hits>1)
df3.agg$five<-(df3.agg$Hic_hits>4)
df3.agg$zero<-(df3.agg$Hic_hits>-1)
df3.agg2 <- aggregate(cbind(zero,one,two,five) ~ Sample+Patient, data=df3.agg,sum)
df3.agg2$oneperc <- df3.agg2$one/df3.agg2$zero
df3.agg2$twoperc <- df3.agg2$two/df3.agg2$zero
df3.agg2$fiveperc <- df3.agg2$five/df3.agg2$zero
write.table(df3.agg2, "Taxa_hic_coverage.txt",row.names=FALSE,quote = F,sep="\t")
df3.agg2$zero <-NULL
df3.agg2$one <- NULL
df3.agg2$two <- NULL
df3.agg2$five <- NULL
colnames(df3.agg2)<-c("Sample","Patient","One+","Two+","Five+")
write.table(stat.desc(df3.agg2),row.names=FALSE,quote = F,"Taxa_hic_coverage_statistics.txt",sep="\t")
df3.agg2.melt <- melt(df3.agg2)
pdf("Taxa_hic_coverage.pdf")
ggplot(df3.agg2.melt, aes(x=Patient, y=value))+
geom_boxplot(outlier.shape = NA)+
geom_point()+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
facet_grid(.~variable)
dev.off() 
