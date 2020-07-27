
library(ggplot2)
library(RColorBrewer)
setwd('/workdir/users/agk85/CDC/newhic/plotting/')
workdir='/workdir/users/agk85/CDC'

args <- commandArgs(trailingOnly = TRUE)

###ggplot theme no background
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())

#load in data
identity=args[1]
df <- read.csv(paste("CDC_percent_capture_", identity,"id.csv",sep=""),row.names=1)
#what's your cutoff? number of links minimum---must have calculated this before---in hic_trans_interactions.R
cutoff=args[2]
df2<- subset(df, grp=="Captured.by.trans.interactions" & cfs==cutoff & type !='trans' & type !='org')
df2$Patient<- as.factor(sapply(strsplit(as.character(df2$sampleids),"\\-"), `[`, 1))
df3 <- subset(df2, Patient == "B357")
df3 <- subset(df2, Patient == "B314")
color_pallete_function <- colorRampPalette(
  colors = c("blue3","burlywood","deepskyblue","darkorange2","aquamarine","darkorchid4","aliceblue"),
  space = "Lab") # Option used when colors do not represent a quantitative scale)

B357_pallete_function <- colorRampPalette(colors=c("darkorange","darkorange4"), space="Lab")
B314_function <- colorRampPalette(colors=c("blue", "blue4"), space="Lab")
num_colors <- length(unique(df2$Patient))
num_B357 <- length(unique(df3$sampleids))
num_B314 <- length(unique(df3$sampleids))
print(num_colors)
colorful_colors<- color_pallete_function(num_colors)
B357_colors <- B357_pallete_function(num_B357)
B314_colors <- B314_function(num_B314)

pdf(paste("Percent_capture_trans_", identity, "id_min_", cutoff, ".pdf",sep=""),height=7, width=10)
ggplot(df3,aes(fill=sampleids,factor(type), mge_percents))+
        geom_bar(position='dodge', stat='identity')+
        geom_text(position=position_dodge(width=.9),cex=2,aes(label=paste(mge_captured, '/\n',mge_total,sep="")), vjust=-0.25)+
        ylim(0,100)+         #max(df2$mge_percents))+
        xlab(label=c("Annotation"))+ 
        ylab(label=c("Percentage of total MGEs"))+
        scale_fill_manual("Timepoints", values =B314_colors)+
        scale_x_discrete(labels=c("arg" = "Antibiotic\nResistance\nGene", "is" = "Insertion\nSequence\nElement","phage"="Phage","plasmid"="Plasmid"))
dev.off()

#from before
colorful_colors <- c("lightgoldenrod","plum","slategray","lemonchiffon3","peru","cyan3","aliceblue")
pdf(paste("Percent_capture_trans_boxplot_", identity, "id_min_", cutoff, ".pdf",sep=""),height=7, width=10,useDingbats=F)
ggplot(df2, aes(type, mge_percents))+
	geom_boxplot(aes(fill=Patient), position="dodge",outlier.shape = NA)+
	geom_point(aes(col=Patient),position = position_jitterdodge(jitter.width = .01, dodge.width = .8))+
	ylim(0,100)+         #max(df2$mge_percents))+
	xlab(label=c("Contig Annotation"))+ 
	ylab(label=c("Percentage of total annotated contigs"))+
	scale_fill_manual("Patients", values = colorful_colors)+
	scale_colour_manual("Patients", values = rep("black", length(df2$sampleids)))+
	scale_x_discrete(labels=c("arg" = "Antibiotic\nResistance\nGene", "is" = "Insertion\nSequence\nElement","phage"="Phage","plasmid"="Plasmid"))
dev.off()


pdf(paste("Percent_capture_trans_boxplot_ungrouped_", identity, "id_min_", cutoff, ".pdf",sep=""),height=7, width=5, useDingbats=F)
ggplot(df2, aes(type, mge_percents))+
        geom_boxplot(position="dodge",outlier.shape = NA)+
        geom_point(aes(col=Patient),position = position_jitter(width=0.2))+
        ylim(0,100)+         #max(df2$mge_percents))+
        xlab(label=c("Contig Annotation"))+
        ylab(label=c("Percentage of total annotated contigs"))+
        #scale_fill_manual("Patients", values = colorful_colors)+
	#theme(legend.position="none")+
        scale_colour_manual("Patients", values=colorful_colors)+ # values = rep("black", length(df2$sampleids)))+
        scale_x_discrete(labels=c("arg" = "Antibiotic\nResistance\nGene", "is" = "Insertion\nSequence\nElement","phage"="Phage","plasmid"="Plasmid"))
dev.off()

