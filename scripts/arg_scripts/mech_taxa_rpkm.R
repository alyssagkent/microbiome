library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(scales)
library(ggplot2)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.grid = element_blank())}
theme_set(theme_nogrid())

#colors
arg.palette <-colorRampPalette(brewer.pal(12,"Set3"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 7
kewl<-col_vector[60:(n+60)]


setwd("/workdir/users/agk85/CDC2/args/abundances")

infile2=read.csv("/workdir/users/agk85/CDC2/args/mapping/bwa_alignments_99_99/arg_v_samp_99_99.txt",header =F,row.names=1)
genes = as.character(unlist(infile2[1,]))
infile3 = infile2[2:nrow(infile2),]
colnames(infile3)<-genes
df3 = t(infile3)
df4<- df3[,order(colnames(df3))]
df4half<-apply(df4, 2, as.numeric)
rownames(df4half) <-rownames(df4)
d<-df4half[rowSums(df4half)<1,]
linkerfile = read.csv("/workdir/users/agk85/CDC2/args/abundances/arg_v_samp_99_99_names_mech_refinder_drug.txt",header=T,sep="\t")
df5<-merge(linkerfile,d,  by.x="Protein",by.y="row.names", all.y=T)
df5$Protein <-NULL
df5$CARD <-NULL
df5$Resfams <-NULL
df5.melt <- melt(df5, id = c("Cluster", "Name", "Mechanism","Sub_mechanism","Drug"))
colnames(df5.melt)<-c("Cluster","Name","Mechanism","Sub_mechanism","Drug","Sample","RPKM")
df5.melt$RPKM<-as.numeric(as.character(df5.melt$RPKM))
df5.melt$Patient <- as.vector(sapply(strsplit(as.character(df5.melt$Sample),"-"), function(x) x[[1]]))



tp_dp <- read.csv("TP_DP.txt",sep="\t",header=T)
tp_dp$Trimmed <-NULL
df6 <- merge(df5.melt, tp_dp,by=c("Sample","Patient"))
df7 <- subset(df6, Sub_mechanism != "Efflux")

drugs=unique(df7$Sub_mechanism)
print(drugs)
for (i in seq(length(drugs))){
drug = drugs[i]
df7pat <- subset(df7,Sub_mechanism==drug)
print(dim(df7pat))
g = ggplot(df7pat,aes(x=DP,y=RPKM,group=Cluster))+
geom_point()+
geom_line(aes(col=Patient))+
facet_grid(Sub_mechanism~Patient,scales="free")+
geom_text(aes(label=ifelse(RPKM>50,as.character(Name),'')),hjust=0,vjust=0)+
scale_color_manual(values=kewl)
name=paste("ARG_abundances_by_mechanism_drug_",drug,".pdf",sep="")
print(name)
pdf(paste("ARG_abundances_by_mechanism_drug_",drug,".pdf",sep=""),useDingbats=F,height=2,width=30)
print(g)
dev.off()
}

patients=unique(df7$Patient)
print(patients)
for (i in seq(length(patients))){
patient = patients[i]
kolor = kewl[i]
df7pat <- subset(df7,Patient==patient)
print(dim(df7pat))
g = ggplot(df7pat,aes(x=DP,y=RPKM,group=Cluster))+
geom_point()+
geom_line(aes(col=Patient))+
facet_grid(Sub_mechanism~Patient,scales="free")+
geom_text(aes(label=ifelse(RPKM>50,as.character(Name),'')),hjust=0,vjust=0)+
scale_color_manual(values=kolor)
name=paste("ARG_abundances_by_mechanism_patient_",patient,".pdf",sep="")
print(name)
pdf(paste("ARG_abundances_by_mechanism_patient_",patient,".pdf",sep=""),useDingbats=F,height=20,width=5)
print(g)
dev.off()
}



g = ggplot(df7,aes(x=DP,y=RPKM,group=Cluster))+
geom_point()+
geom_line(aes(col=Patient))+
facet_grid(Sub_mechanism~Patient,scales="free")+
geom_text(aes(label=ifelse(RPKM>50,as.character(Name),'')),hjust=0,vjust=0)+
scale_color_manual(values=kewl)
pdf(paste("ARG_abundances_by_mechanism.pdf",sep=""),useDingbats=F,height=18,width=30)
print(g)
dev.off()


df8 <- subset(df6, Drug %in% c("piperacillin-tazobactam","levofloxacin","sulfamethoxazole","trimethoprim"))
pdf("ARG_abundances_by_drug.pdf",height=10,width=20,useDingbats=F)
ggplot(df8,aes(x=DP,y=RPKM,group=Cluster))+
geom_point()+
geom_line(aes(col=Patient))+
facet_grid(Drug~Patient,scales="free")+
geom_text(aes(label=ifelse(RPKM>0.5,as.character(Name),'')),hjust=0,vjust=0)+
scale_color_manual(values=kewl)
dev.off()

pdf("ARG_abundances_by_drug_nonames.pdf",height=10,width=20,useDingbats=F)
ggplot(df8,aes(x=DP,y=RPKM,group=Cluster))+
geom_point()+
geom_line(aes(col=Patient))+
facet_grid(Drug~Patient,scales="free")+
scale_color_manual(values=kewl)
dev.off()

