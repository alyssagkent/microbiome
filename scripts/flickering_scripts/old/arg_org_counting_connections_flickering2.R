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

folder = 'CDC2'
infolder = '/workdir/users/agk85/CDC2/bins/flickering2'
inhandle = '/workdir/users/agk85/CDC2/bins/flickering2/arg_org_hic_separated_all_5.tbl'
genetype = 'arg'
minreads = '5'
level = 'all'

args <- commandArgs(trailingOnly = TRUE)
infolder = args[1]
inhandle = args[2]
genetype = args[3]
minreads = args[4]
level = args[5]


setwd(infolder)
indata = read.csv(inhandle, header=T, sep="\t",row.names=1)

colordf = indata[,grep("color",colnames(indata))]

p1 <-as.vector(sapply(strsplit(colnames(colordf),"color"), function(x) x[[1]]))
patients <-as.vector(sapply(strsplit(p1, "[.]"), function(x) x[[1]]))
pats<- unique(patients)
colordf$Cluster<-sub2$Cluster
colordf$Organism<-sub2$Organism
colordf$ARG_name<-sub2$ARG_name
#4,7,4,4,3,6,5,5,5
starts <- c(1,5,12,16,20,23,29,34,39)
stops <- c(4, 11, 15,19,22,28,33,38,43)
avg_metric<-list()
pass_metric<-list()
all_metric<-list()
for (i in seq(length(pats))){
	patient<- pats[i]
	print(patient)
	start <-starts[i]
	stop <- stops[i]
	poi<- colordf[,start:stop]
fours <-rowSums(as.matrix((poi > 3) + 0))
threes <- rowSums(as.matrix((poi  > 2) + 0))-fours
poi$fours <-fours
poi$threes<-threes
poi$Cluster <- colordf$Cluster
poi$ARG_name <- colordf$ARG_name
poi$Organism <- colordf$Organism
poiconnect <- subset(poi, poi$fours>0)
poiconnect$proportion <- poiconnect$fours /(poiconnect$fours + poiconnect$threes)
#get which ones are 1 and which are less than one (can't be more than .99 because only max 7 timepoints)
poiconnect$passfail<-(poiconnect$proportion>.99) + 0
print(mean(poiconnect$proportion))
print(sum(poiconnect$passfail==1)/length(poiconnect$passfail))
print(sum(poiconnect$fours)/(sum(poiconnect$fours)+sum(poiconnect$threes)))
write.csv(poiconnect, paste(patient, "_poiconnect.csv",sep=""))
avg_metric[[i]]<-mean(poiconnect$proportion)
pass_metric[[i]]<-sum(poiconnect$passfail==1)/length(poiconnect$passfail)
all_metric[[i]]<-sum(poiconnect$fours)/(sum(poiconnect$fours)+sum(poiconnect$threes))
print(sum(poiconnect$fours))
print(sum(poiconnect$threes))

}

#####
d<- data.frame(pats, unlist(avg_metric), unlist(pass_metric),unlist(all_metric))
colnames(d)<-c("Patients","Average","Pass","Alltogether")

pdf(paste("Flickering_",genetype,"_org_",minreads, "_", level,"_patient_passing_base.pdf",sep=""))
ggplot(data=d,aes(Patients,Pass))+
	geom_bar(stat="identity")+
	ylim(0,1)
dev.off()

pdf(paste("Flickering_",genetype,"_org_",minreads, "_", level,"_patient_alltogether_base.pdf",sep=""))
ggplot(data=d,aes(Patients,Alltogether))+
	geom_bar(stat="identity")+
	ylim(0,1)
dev.off()

pdf(paste("Flickering_",genetype,"_org_",minreads, "_", level,"_patient_average_base.pdf",sep=""))
ggplot(data=d,aes(Patients,Average))+
	geom_bar(stat="identity")+
	ylim(0,1)
dev.off()

