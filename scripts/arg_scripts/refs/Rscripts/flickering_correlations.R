#arg_org_counting_connections_vs_not.R
library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(scales)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}

theme_set(theme_nogrid())

##################################################
#color generation
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 7
kewl<-col_vector[1:(n+1)]

setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering")
indata = read.csv("arg_org_hic_separated_95_98_1_2.tbl", header=T, sep="\t",row.names=1)
#indata = read.csv("arg_org_hic_separated_95_98_0_2.tbl", header=T, sep="\t",row.names=1)

sub2<-indata
orgdf = sub2[,grep("_org",colnames(sub2))]
argdf = sub2[,grep("_arg",colnames(sub2))]
connectdf = sub2[,grep("_connect",colnames(sub2))]
colordf = sub2[,grep("color",colnames(sub2))]

cd<-colordf[colordf>3]

test <-head(colordf)
#this worked ont he cluster but not on my computer idk-why

colnames(colordf)
patients <-as.vector(sapply(strsplit(colnames(colordf),"[.]"), function(x) x[[1]]))
pats<- unique(patients)
colordf$Cluster<-sub2$Cluster
colordf$Organism<-sub2$Organism
colordf$ARG_name<-sub2$ARG_name

starts <- c(1,5,12,17,21,24,30)
stops <- c(4, 11, 16,20,23,29,32)
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
}

#####
d<- data.frame(pats, unlist(avg_metric), unlist(pass_metric),unlist(all_metric))
colnames(d)<-c("Patient","Average","Pass","Alltogether")

#pdf("Flickering_argorg_patient_alltogether_base.pdf")
ggplot(data=d,aes(Patient,Alltogether))+
	geom_bar(stat="identity")+
	ylim(0,1)
#dev.off()

#sample	patient mgm_reads mgm_assembly hic_reads hic_trans_reads

#then plot flickering correlation figures
sub_d <- subset(d, select=c(Patient, Alltogether))
metadata <- read.csv('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/metacounts_combined.txt',header=T,sep="\t")
combo_d <- merge(sub_d, metadata, by="Patient", all.y=T)

cl = data.frame(combo_d$Patient, kewl)
d_summary <- ddply(combo_d,~Patient,summarise,
Mgm_reads_avg=mean(Mgm_reads),Mgm_reads_min=min(Mgm_reads), Mgm_reads_total=sum(Mgm_reads),
Hic_reads_avg=mean(Hic_reads),Hic_reads_min=min(Hic_reads), Hic_reads_total=sum(Hic_reads),
Mgm_assembly_avg=mean(Mgm_assembly),Mgm_assembly_min=min(Mgm_assembly), Mgm_assembly_total=sum(Mgm_assembly),
Hic_trans_reads_avg=mean(Hic_trans_reads),Hic_trans_reads_min=min(Hic_trans_reads), Hic_trans_reads_total=sum(Hic_trans_reads))
d_sum <- merge(sub_d, d_summary, by="Patient")
d_sum$timepoints <- c(4,7,5,4,3,6,3)
cool <- kewl
names(cool)<-d_sum$Patient
a = ggplot(d_sum, aes(x=Mgm_reads_avg, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Mgm_reads_avg, d_sum$Alltogether),2), x = 35000000 , y = .8 , size = 3)
b = ggplot(d_sum, aes(x=Mgm_reads_min, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Mgm_reads_min, d_sum$Alltogether),2), x = 18000000, y = .8 , size = 3)
c = ggplot(d_sum, aes(x=Mgm_reads_total, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Mgm_reads_total, d_sum$Alltogether),2), x = 150000000, y = .8 , size = 3)

d = ggplot(d_sum, aes(x=Hic_reads_avg, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Hic_reads_avg, d_sum$Alltogether),2), x = 15000000, y = .8 , size = 3)
e = ggplot(d_sum, aes(x=Hic_reads_min, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Hic_reads_min, d_sum$Alltogether),2), x = 5000000, y = .8 , size = 3)
f = ggplot(d_sum, aes(x=Hic_reads_total, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Hic_reads_total, d_sum$Alltogether),2), x = 60000000, y = .8 , size = 3)

g = ggplot(d_sum, aes(x=Mgm_assembly_avg, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Mgm_assembly_avg, d_sum$Alltogether),2), x = 100000000, y = .8 , size = 3)
h = ggplot(d_sum, aes(x=Mgm_assembly_min, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Mgm_assembly_min, d_sum$Alltogether),2), x = 70000000, y = .8 , size = 3)
i = ggplot(d_sum, aes(x=Mgm_assembly_total, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Mgm_assembly_total, d_sum$Alltogether),2), x = 1300000000, y = .8 , size = 3)

j = ggplot(d_sum, aes(x=Hic_trans_reads_avg, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Hic_trans_reads_avg, d_sum$Alltogether),2), x = 25000, y = .8 , size = 3)
k = ggplot(d_sum, aes(x=Hic_trans_reads_min, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Hic_trans_reads_min, d_sum$Alltogether),2), x = 5000, y = .8 , size = 3)
l = ggplot(d_sum, aes(x=Hic_trans_reads_total, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none",axis.text.x=element_text(size=4))+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$Hic_trans_reads_total, d_sum$Alltogether),2), x = 100000, y = .8 , size = 3)
m = ggplot(d_sum, aes(x=timepoints, y=Alltogether, col=Patient))+ scale_color_manual(values=cool)+ geom_point()+
theme(legend.position="none")+scale_x_continuous(labels = comma)+
annotate("text", label = round(cor(d_sum$timepoints, d_sum$Alltogether),2), x = 6.5, y = .8 , size = 3)
pdf("Flickering_correlations.pdf")
grid.arrange(a,b,c,d,e,f,g,h,i,j,k,l,m, ncol=3, nrow=5, widths=c(6,6,6),heights=c(6,6,6,6,6))
dev.off()


