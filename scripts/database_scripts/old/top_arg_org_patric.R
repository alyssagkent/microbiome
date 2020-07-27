library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library("devtools")
library(igraph)
set.seed(42)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x, palette){
  cols <- colorRamp(palette(100))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}
##################################################
#color generation
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n<-14
mech.col<-col_vector[59:(n+59)]
#pie(rep(1,n), col=mech.col)
n <- 7
kewl<-col_vector[1:(n+1)]
#pie(rep(1,n), col=kewl)
n<-32
family.col<-sample(col_vector, n,replace = FALSE)
#pie(rep(1,n), col=family.col)
n<-7
phylum.col<-sample(col_vector, n,replace = FALSE)
phylum.col<-col_vector[47:(n+47)]
#pie(rep(1,n), col=phylum.col)

###################################################################
setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/taxonomic_distributions")
###################################################################
library(reshape2)
#import arg rpkm
infile2=read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95.txt",header =F,row.names=1)
genes = as.character(unlist(infile2[1,]))
infile3 = infile2[2:nrow(infile2),]
colnames(infile3)<-genes
df3 = t(infile3)
df4<- df3[,order(colnames(df3))]
linkerfile = read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_v_samp_95_95_names_mech.txt",header=T,sep="\t")
df5<-merge(linkerfile,df4,  by.x="Protein",by.y="row.names", all.y=T)
df5$Protein <-NULL
df5$CARD <-NULL
df5$Resfams <-NULL
df5.melt <- melt(df5, id = c("Cluster", "Name", "Mechanism","Sub_mechanism"))
colnames(df5.melt)<-c("Cluster","Name","Mechanism","Sub_mechanism","Sample","RPKM")
df5.melt$RPKM<-as.numeric(as.character(df5.melt$RPKM))
df5.melt$Patient <- as.vector(sapply(strsplit(as.character(df5.melt$Sample),"-"), function(x) x[[1]]))
#aggregate down to the Patient just in case
df5.agg <- aggregate(RPKM ~ Cluster + Name + Mechanism + Sub_mechanism + Patient, data = df5.melt, max)
####################################################################################################################

depth = 1
#import the cluster info
arginfo= read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_v_samp_95_95_names_mech.txt", sep="\t",header = T)
indata = read.csv(paste("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_org_hic_cliques_95_98_",depth, "_2.tbl",sep=""), header=T, sep="\t",row.names=1)

indata_file = read.csv(paste("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_org_hic_cliques_95_98_",depth, "_2.tbl",sep=""), header=F, sep="\t",row.names=1)
taxa = as.character(unlist(indata_file[1,]))

sub2<- indata[order(indata$ARG_name, indata$Cluster),] #indata$Sample, do first if want to reorder
#sub2 <- subset(sub2, Top_ARG!='NA') 
sub4 <- sub2
connectdf = sub4[5:ncol(sub4)]
connectsum <- rowSums(connectdf)
#everything will either be normal or multitaxa
sub5 <- subset(sub4,connectsum>0)
connectdf <-sub5[5:ncol(sub5)]
connectdf <- connectdf[,colSums(connectdf)>0]
n=3
levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))

subdata = data.matrix(connectdf)
colors <- c("white","black")
heat.col = colors[0:max(subdata)+1]
n <- max(subdata)
namelist = c("not connected","connected")
names <- namelist[0:max(subdata)+1]
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}
sub5$patients <- as.vector(sapply(strsplit(as.character(sub5$Sample),"-"), function(x) x[[1]]))
############################
patients <- sub5$patients
pats = unique(patients)[order(unique(patients))]
cols <- colorRampPalette(kewl)
patient.col = cols(length(pats))

###################################################
#top ARGs
###################################################
#another follow-up question
#where are the top-10-ARGs who are they residing in??
toparg_clusters <-read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/topargs/all_clusters_topargs_mergednames.txt",header=T,sep="\t")
ta_clusters <- unlist(toparg_clusters)
top10<-merge(sub5,toparg_clusters,by="Cluster")
topgeneral<-read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/top_args_categories.csv", header=T)
known_cat<-read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/topargs/known_categories.txt",sep="\t",header=F)
colnames(known_cat)<-c("Mergedname","Category","Trim_name")
top10 <-merge(known_cat,top10,by="Mergedname")

sub_argcount <- top10
sub_argcount$connectsum <-NULL
sub_argcount$Trim_name<-NULL
sub_argcount$Sample <-NULL
sub_argcount$Confidence <-NULL
sub_argcount$ARG_name<-NULL
sub_argcount$Specific_name<-NULL

df_melt2 <- melt(sub_argcount, id = c("patients", "Cluster", "Top_ARG","Mergedname","Category"))
df_melt.pos <-subset(df_melt2, value>0)
df_melt.pos$levels <-as.vector(sapply(strsplit(as.character(df_melt.pos$variable),"[.][.]"), function(x) list(trimws(x))))
df_melt.pos$family <- sapply(df_melt.pos$levels, function(x) paste(x[1:5], collapse=";"))
df_melt.pos$genus <- sapply(df_melt.pos$levels, function(x) paste(x[1:6], collapse=";"))
df_melt.pos$species <- sapply(df_melt.pos$levels, function(x) paste(x[1:7], collapse=";"))
df_melt.pos$familyname <- sapply(df_melt.pos$levels, function(x) paste(x[5:5], collapse=";"))
df_melt.pos$genusname <- sapply(df_melt.pos$levels, function(x) paste(x[6:6], collapse=";"))
df_melt.pos$speciesname <- sapply(df_melt.pos$levels, function(x) paste(x[7:7], collapse=";"))

#family
df.arg.agg.family <-aggregate(value ~patients+Cluster+Mergedname+Category+family+familyname, data =  df_melt.pos, sum)
df.arg.agg.family$value <- (df.arg.agg.family$value>0) + 0
colnames(df.arg.agg.family)<-c("Patient","Cluster","Mergedname","Category","Family","Familyname","Count")
df.arg.agg.family.reduced<-aggregate(Count~Patient+Family + Familyname,data=df.arg.agg.family,sum)
colnames(df.arg.agg.family.reduced)<-c("Patient","Family","Familyname","Count")

#this will collapse on the patients and keep the clusters unique
df.arg.agg.family.patient.reduced<-aggregate(Count~Family+Familyname+Cluster+Mergedname+Category,data=df.arg.agg.family,sum)
df.arg.agg.family.patient.reduced$Count <- (df.arg.agg.family.patient.reduced$Count>0) + 0
colnames(df.arg.agg.family.patient.reduced)<-c("Family","Familyname","Cluster","Mergedname","Category","Count")
pdf(paste("Unique_TopARGs_categories_family_clusterpatientcollapse.pdf",sep=""),height=20, width=15,useDingbats=F)
plot19 <-ggplot(df.arg.agg.family.patient.reduced, aes(x = Family, y = Count, fill = Category))+
	geom_bar(stat="identity")+
	xlab("Taxa")+
	ylab("ARG Count")+
	#scale_fill_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot19)
dev.off()

patric <- read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/taxonomic_distributions/Patric_toparg_taxonomies.txt",header=F, sep="\t")
colnames(patric) <- c("Cluster","Patric_taxonomy")

toppatric <- merge(patric, toparg_clusters, by="Cluster")
toppatric <- merge(known_cat, toppatric, by="Mergedname")
toppatric$Trim_name <- NULL

df_patric <- toppatric
df_patric$levels <-as.vector(sapply(strsplit(as.character(df_patric$Patric_taxonomy),"; "), function(x) list(trimws(x))))
df_patric$Family <- sapply(df_patric$levels, function(x) paste(x[1:5], collapse=";"))
df_patric$Familyname <- sapply(df_patric$levels, function(x) paste(x[5:5], collapse=";"))
df_patric$Count <- rep(1,nrow(df_patric))
df_patric$levels <-NULL
df_patric$Patric_taxonomy <-NULL
#family
#df.arg.agg.family <-aggregate(value ~Cluster+Mergedname+Category+family+familyname, data =  df_patric, sum)
#df.arg.agg.family$value <- (df.arg.agg.family$value>0) + 0
#colnames(df.arg.agg.family)<-c("Cluster","Mergedname","Category","Family","Familyname","Count")
#df.arg.agg.family.reduced<-aggregate(Count~Family + Familyname,data=df.arg.agg.family,sum)
#colnames(df.arg.agg.family.reduced)<-c("Family","Familyname","Count")

#this will collapse on the patients and keep the clusters unique
df.arg.agg.family.patient.reduced<-aggregate(Count~Family+Familyname+Cluster+Mergedname+Category,data=df_patric,sum)
df.arg.agg.family.patient.reduced$Count <- (df.arg.agg.family.patient.reduced$Count>0) + 0
colnames(df.arg.agg.family.patient.reduced)<-c("Family","Familyname","Cluster","Mergedname","Category","Count")
pdf(paste("Unique_TopARGs_Patric_categories_family_clustercollapse.pdf",sep=""),height=20, width=15,useDingbats=F)
plot20 <-ggplot(df.arg.agg.family.patient.reduced, aes(x = Family, y = Count, fill = Category))+
        geom_bar(stat="identity")+
        xlab("Taxa")+
        ylab("ARG Count")+
        #scale_fill_manual(label=patients, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot20)
dev.off()


