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

setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering")
indata = read.csv("arg_org_hic_separated_95_98_1_2.tbl", header=T, sep="\t",row.names=1)
indata = read.csv("arg_org_hic_separated_95_98_0_2.tbl", header=T, sep="\t",row.names=1)

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
colnames(d)<-c("Patients","Average","Pass","Alltogether")

pdf("Flickering_argorg_patient_passing_base.pdf")
ggplot(data=d,aes(Patients,Pass))+
	geom_bar(stat="identity")+
	ylim(0,1)
dev.off()

pdf("Flickering_argorg_patient_alltogether_base.pdf")
ggplot(data=d,aes(Patients,Alltogether))+
	geom_bar(stat="identity")+
	ylim(0,1)
dev.off()

pdf("Flickering_argorg_patient_average_base.pdf")
ggplot(data=d,aes(Patients,Average))+
	geom_bar(stat="identity")+
	ylim(0,1)
dev.off()

sub2$orgsum <- rowSums(orgdf)
sub2$argsum <- rowSums(argdf)
sub2$connectsum <- rowSums(connectdf)

sub3 <- subset(sub2, orgsum>1 )
sub4 <- subset(sub3, argsum>0 & ARG_name=='CMY-LAT-MOX-ACT-MIR-FOX')
n=5
levels <-as.vector(sapply(strsplit(as.character(sub4$Organism),"; "), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
sub4$Taxa <- groupvar
sub5 <- subset(sub4, Taxa=='f__Enterobacteriaceae')

colordf = sub5[,grep("color",colnames(sub5))]
subdata = data.matrix(colordf)
colors <- c("black","white","red","blue","orange")
heat.col = colors[0:max(subdata)+1]
n <- max(subdata)
namelist = c("neither","org","arg","not connected","connected")
names <- namelist[0:max(subdata)+1]
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="single")}


#pdf(paste("heatmap_arg_org_depth1_top10_none.pdf",sep=""), height=10, width=15,useDingbats=FALSE)
heatmap.2(subdata, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="none",	Rowv=F,Colv=F,key=F,symbreaks=FALSE, margins=c(10,15),
	symkey=F, density.info="none", trace="none", labCol=colnames(subdata),
	labRow=sub5$Organism,	cexRow=.4, cexCol=.8, col=heat.col)
 legend("left",legend=names,
	fill=c(heat.col), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

#dev.off()


#pdf(paste("heatmap_arg_org_depth1_top10_one.pdf",sep=""), height=10, width=15,useDingbats=FALSE)
heatmap.2(subdata, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="row",	Rowv=T,Colv=F,key=F,symbreaks=FALSE, margins=c(10,15),
	symkey=F, density.info="none", trace="none", labCol=colnames(subdata),
	labRow=sub5$Organism,	cexRow=.4, cexCol=.8, col=heat.col)
 legend("left",legend=names,
	fill=c(heat.col), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
#dev.off()

 
 
########################get samps vs. args
 
argdf = sub2[,grep("_arg",colnames(sub2))]
sub_argcount <- argdf
sub_argcount$Top_ARG <-NULL
sub_argcount$Organism <-NULL
df_1 <- cbind(sub2$Cluster, sub2$ARG_name, sub_argcount)
samps<-unlist(as.vector(sapply(strsplit(colnames(sub_argcount),"_"), function(x) list(trimws(x)[[1]]))))
pats <-unlist(as.vector(sapply(strsplit(colnames(sub_argcount),"[.]"), function(x) list(trimws(x)[[1]]))))

colnames(df_1)<- c("Cluster", "ARG_name",samps)
df_melt2 <- melt(df_1, id = c("Cluster", "ARG_name"))
colnames(df_melt2)
df_melt.pos <-subset(df_melt2, value>0)
df.arg.agg <-aggregate(value ~ variable+Cluster+ARG_name, data =  df_melt.pos, sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg$value <- (df.arg.agg$value>0)+0
colnames(df.arg.agg)<-c("Sample","Cluster","ARG_name","Presence")

df.arg.agg$Patient <- unlist(as.vector(sapply(strsplit(as.character(df.arg.agg$Sample),"[.]"), function(x) list(trimws(x)[[1]]))))
df.arg.patient <-aggregate(Presence ~ Patient+Cluster+ARG_name, data =  df.arg.agg, sum)
df.arg.patient$Presence <- (df.arg.patient$Presence>0)+0
df.arg.shared<-aggregate(Presence ~ Cluster+ ARG_name, data =  df.arg.patient, sum)
#what percentage of ARG clusters have more than one patient?
sum(df.arg.shared$Presence>1)/nrow(df.arg.shared)

df.arg.agg.temp<-df.arg.agg
df.arg.agg.temp$Patient<-NULL
arg_sample <- dcast(df.arg.agg.temp, Cluster + ARG_name ~ Sample, sum)
rownames(arg_sample)<-arg_sample$Cluster
arg_sample$Cluster<-NULL

###plot that figure
arg.colors <- colorRampPalette(brewer.pal(12,"Set3"))
arg.col = arg.colors(length(unique(arg_sample$ARG_name)))
rsc=arg.col[factor(arg_sample$ARG_name)]
args <- arg_sample$ARG_name
arg_sample$ARG_name<-NULL
pat.colors <- colorRampPalette(brewer.pal(9,"Set1"))
patients <-as.vector(sapply(strsplit(colnames(arg_sample),"[.]"), function(x) x[[1]]))
pat.col = pat.colors(length(unique(patients)))
csc=pat.col[factor(patients)]

a<-unique(data.frame(args, rsc))
d<-unique(data.frame(patients,csc))


mat<-arg_sample
#mat$ARG_name<-NULL
mat2<-as.matrix(mat)
my_palette <- c("white", "red")
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}

pdf(paste("heatmap_arg_samp_presence_none.pdf",sep=""), height=10, width=15,useDingbats=FALSE)
heatmap.2(mat2, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="none",	Rowv=F,Colv=T,key=F,symbreaks=FALSE, margins=c(10,18),
	symkey=F, density.info="none", trace="none", labCol=colnames(mat2),
	labRow=NA,	cexRow=.4, cexCol=1.5, col=my_palette, ColSideColors=csc, RowSideColors=rsc)

legend("topleft",legend=c(as.character(d$patients),"absence","presence"),
	fill=c(as.character(d$csc),"white", "red"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

#legend("topleft",legend=c( as.character(a$newnames), as.character(d$levels),"0",as.character(max(mat))),
#	fill=c(as.character(a$rsc), as.character(d$csc),"white", "red"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

#############aggregated over patients instead
df.arg.agg.temp<-df.arg.patient
arg_patient <- dcast(df.arg.agg.temp, Cluster + ARG_name ~ Patient, sum)
rownames(arg_patient)<-arg_sample$Cluster
arg_patient$Cluster<-NULL

###plot that figure
arg.colors <- colorRampPalette(brewer.pal(12,"Set3"))
arg.col = arg.colors(length(unique(arg_patient$ARG_name)))
rsc=arg.col[factor(arg_patient$ARG_name)]
args <- arg_patient$ARG_name
arg_patient$ARG_name<-NULL
pat.colors <- colorRampPalette(brewer.pal(9,"Set1"))
patients <-colnames(arg_patient)   #as.vector(sapply(strsplit(colnames(arg_patient),"[.]"), function(x) x[[1]]))
pat.col = pat.colors(length(unique(patients)))
csc=pat.col[factor(patients)]

a<-unique(data.frame(args, rsc))
d<-unique(data.frame(patients,csc))


mat<-arg_patient
#mat$ARG_name<-NULL
mat2<-as.matrix(mat)
my_palette <- c("white", "red")
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}



pdf(paste("heatmap_arg_samp_presence_patient_both.pdf",sep=""), height=10, width=15,useDingbats=FALSE)
heatmap.2(mat2, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="both",	Rowv=T,Colv=T,key=F,symbreaks=FALSE, margins=c(10,18),
	symkey=F, density.info="none", trace="none", labCol=colnames(mat2),
	labRow=NA,	cexRow=.4, cexCol=1.5, col=my_palette, ColSideColors=csc, RowSideColors=rsc)

legend("topleft",legend=c(as.character(d$patients),"absence","presence"),
	fill=c(as.character(d$csc),"white", "red"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()


