#arg_org_counting_connections_vs_not.R
library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())



setwd('/workdir/users/agk85/CDC2/bins/flickering')

indata = read.csv('flickering.txt', header=T, sep="\t")
patient<- c("B314","B316","B320","B331","B335","B357","B370","US3","US8")
numsamps<-c(4,7,4,4,3,6,5,5,5)
num_samples = data.frame(patient,numsamps)
df <- merge(indata, num_samples, by="patient")
df$patn <- paste(df$patient, " (n=",df$numsamps,")",sep="")

pdf(paste("Flickering_alltogether.pdf",sep=""),width=7,height=6)
ggplot(data=df,aes(patn,consistency))+
	geom_bar(stat="identity")+
	ylim(0,1)+
	xlab("Patients")+
	theme(axis.text.x = element_text(angle = 90, hjust = 1))+
	facet_grid(genetype~minthresh)
dev.off()

