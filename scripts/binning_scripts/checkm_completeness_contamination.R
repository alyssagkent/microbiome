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
    theme(panel.grid = element_blank())}
theme_set(theme_nogrid())

setwd("/workdir/users/agk85/CDC2/das")

arg.palette <-colorRampPalette(brewer.pal(12,"Set3"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 9
kewl<-col_vector[60:(n+60)]

header = c("Bin","Taxa","Completeness","Contamination","Hetero","Size")
infile = read.csv("all.stats", sep="\t",header=F)
colnames(infile)<-header
df <- infile
#infile$Sample = as.vector(sapply(strsplit(as.character(infile$Bin), "[.]"),function(x) x[[1]]))
pdf("DAS_completeness_contamination_length.pdf",useDingbats=FALSE)
ggplot(df, aes(x=Completeness,y=Contamination))+
geom_point()+
geom_hline(yintercept=10, linetype="dashed", color="black", size=1)
dev.off()

#ordered plot

dfordered = df[order(-df$Completeness),]
dfordered$Rank = seq(1,nrow(dfordered))
pdf("DAS_completeness_contamination_ordered.pdf",useDingbats=FALSE)
ggplot(dfordered, aes(x=Rank))+
geom_point(aes(y=Completeness),col="red")+
geom_point(aes(y=Contamination),col="blue")+
theme( axis.line.y.right = element_line(color = "blue"), 
       axis.ticks.y.right = element_line(color = "blue"),
	axis.text.y.right = element_text(color = "blue"))+
theme( axis.line.y.left = element_line(color = "red"),
       axis.ticks.y.left = element_line(color = "red"),
	axis.text.y.left = element_text(color = "red"))+
geom_hline(yintercept=10, linetype="dashed", color="black", size=1)+
scale_y_continuous(name="Completeness",sec.axis = sec_axis(~ .,name = "Contamination"))
dev.off()








df2 <- subset(df, Contamination<10)
dim(df2)

df2$Patient = sapply(strsplit(as.character(df2$Bin),"-"), function(x) x[[1]])
df2$Taxa <- NULL
df2$Hetero<-NULL
df2.melt<-melt(df2, by=c(Bin, Patient))

pdf("DAS_metrics_by_patient.pdf",height=8,width=10,useDingbats=FALSE)
ggplot(df2.melt, aes(x=Patient,y=value))+
geom_violin(aes(fill=Patient))+
facet_wrap(~variable,nrow=3,scales="free_y",strip.position = "left",
labeller = as_labeller(c(Completeness="Completeness (%)", Contamination = "Contamination (%)", Size = "Bin length (bp)")))+
labs(x="Patients",y=NULL)+
scale_fill_manual(values=kewl)+
scale_y_continuous(labels = comma)
dev.off()

pdf("DAS_metrics_by_patient_boxplot_quality.pdf",height=8,width=10,useDingbats=FALSE)
ggplot(df2.melt, aes(x=Patient,y=value))+
geom_boxplot(aes(fill=Patient))+
geom_jitter(shape=16, position=position_jitter(0.2))+
#geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
facet_wrap(~variable,nrow=3,scales="free_y",strip.position = "left",
labeller = as_labeller(c(Completeness="Completeness (%)", Contamination = "Contamination (%)", Size = "Bin length (bp)")))+
labs(x="Patients",y=NULL)+
scale_fill_manual(values=kewl)+
scale_y_continuous(labels = comma)
dev.off()


library(psych)
hq <- subset(df2, Contamination<5 & Completeness>90)
mq <- subset(df2, Contamination<10 & Completeness>=50)
lq <- subset(df2, Contamination<10 & Completeness<50)

b = describeBy(hq, "Patient")
c = describeBy(mq, "Patient")
d = describeBy(lq, "Patient")

hqstats = rbind(b$B314,b$B316,b$B320,b$B331,b$B335,b$B357,b$B370,b$US3,b$US8)
mqstats = rbind(c$B314,c$B316,c$B320,c$B331,c$B335,c$B357,c$B370,c$US3,c$US8)
lqstats = rbind(d$B314,d$B316,d$B320,d$B331,d$B335,d$B357,d$B370,d$US3,d$US8)

write.csv(hqstats, "HQ_bins_stats.txt")
write.csv(mqstats, "MQ_bins_stats.txt")
write.csv(lqstats, "LQ_bins_stats.txt")
