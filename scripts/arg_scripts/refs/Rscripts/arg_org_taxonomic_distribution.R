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

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
##################################################
#color generation
library(RColorBrewer)
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
setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/enteros")
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
sub6 <- sub5
sub6$Sample<-NULL
sub6$Top_ARG<-NULL
sub6$ARG_name<-as.character(sub6$ARG_name)
############################
#counting number of ARGs that are attached to an entero vs. not
patients <- sub6$patients
pats = unique(patients)[order(unique(patients))]
cols <- colorRampPalette(kewl)
patient.col = cols(length(pats))

#arg count
sub_argcount <- sub6
sub_argcount$connectsum <- NULL
df_melt2 <- melt(sub_argcount, id = c("patients", "Cluster", "ARG_name"))
df_melt.pos <-subset(df_melt2, value>0)
df_melt.pos$levels <-as.vector(sapply(strsplit(as.character(df_melt.pos$variable),"[.][.]"), function(x) list(trimws(x))))
df_melt.pos$family <- sapply(df_melt.pos$levels, function(x) paste(x[5:5], collapse=";"))
df_melt.pos$familyname <- sapply(df_melt.pos$levels, function(x) paste(x[1:5], collapse=";"))
#get rid of things that are above family
df_melt.below<-subset(df_melt.pos, family!="f__")


#aggregate keeping clusters, patients, and taxa
df.arg.agg <-aggregate(value ~ patients +  Cluster+variable + family + familyname, data=df_melt.below, FUN=sum)
df.arg.agg$value <- (df.arg.agg$value>0)+0
colnames(df.arg.agg)<-c("patients","cluster","variable","family","familyname","value")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.fam <-aggregate(value ~ cluster + patients +family+familyname, data =  df.arg.agg, sum)
colnames(df.arg.agg.fam)<-c("cluster","patients","family","familyname","value")
df.arg.agg.fam$value <- (df.arg.agg.fam$value>0)+0
df.arg.agg.fam.reduced <-aggregate(value ~ patients +family+familyname, data =  df.arg.agg.fam, sum)
df.arg.agg.fam.reduced$entero <- ifelse(as.character(df.arg.agg.fam.reduced$family)== "f__Enterobacteriaceae","Enterobacteriaceae","Not Enterobacteriaceae")

pdf(paste("Unique_ARGs_family_per_patient_comparison_allargs_clustercollapse.pdf",sep=""),height=6, width=4,useDingbats=F)
g <-ggplot(df.arg.agg.fam.reduced, aes(factor(df.arg.agg.fam.reduced$entero),df.arg.agg.fam.reduced$value))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(g)
dev.off()
a <- aov(value~entero, data =df.arg.agg.fam.reduced)

###############################################################################
#taking this df5.agg (which is aggregated down to the patient by maxing the abundances)
#how about looking at arg-rpkm vs. enterobacteriaceae for the ones that are linked
#df.arg.agg.fam.merge <- merge(df.arg.agg.fam, df5.agg, by.x=c("cluster", "patients"), by.y=c("Cluster","Patient"), all.x=T)
#df.arg.agg.fam.merge$Name<-NULL
#df.arg.agg.fam.merge$Mechanism<-NULL
#df.arg.agg.fam.merge$Sub_mechanism<-NULL
#df.arg.agg.fam.merge.reduced <-aggregate(. ~ patients +family+familyname, data =  df.arg.agg.fam.merge, sum)
#df.arg.agg.fam.merge.reduced$entero <- ifelse(as.character(df.arg.agg.fam.merge.reduced$family)== "f__Enterobacteriaceae","Enterobacteriaceae","Not Enterobacteriaceae")
#pdf("Unique_ARGs_family_patient_comparison_allargs_RPKM_clustercollapse.pdf",useDingbats=F, height=8, width=4)
#g <-ggplot(df.arg.agg.fam.merge.reduced, aes(factor(df.arg.agg.fam.merge.reduced$entero),df.arg.agg.fam.merge.reduced$RPKM))+
#	geom_boxplot(outlier.shape=NA)+
#	geom_point(aes(col=patients), position=position_jitter())+
#	xlab("Taxa")+
	ylab("ARG RPKM summed")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(g)
dev.off()




##################
pdf(paste("Unique_ARGs_family_per_patient_allargs_clustercollapse.pdf",sep=""),height=10, width=15,useDingbats=F)
h <-ggplot(df.arg.agg.fam.reduced, aes(factor(df.arg.agg.fam.reduced$familyname),df.arg.agg.fam.reduced$value))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(h)
dev.off()

##################################################ARG abundance#############################################################
#pdf(paste("Unique_ARGs_family_per_patient_allargs_RPKM_clustercollapse.pdf",sep=""),height=10, width=15,useDingbats=F)
#h.rpkm <-ggplot(df.arg.agg.fam.merge.reduced, aes(factor(df.arg.agg.fam.merge.reduced$familyname),df.arg.agg.fam.merge.reduced$RPKM))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG RPKM summed")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(h.rpkm)
dev.off()


#follow up plot
#which bacteria make up the entero arg connections
sub_argcount <- sub6
sub_argcount$connectsum <- NULL
df_melt2 <- melt(sub_argcount, id = c("patients", "Cluster", "ARG_name"))
df_melt.pos <-subset(df_melt2, value>0)
df_melt.pos$levels <-as.vector(sapply(strsplit(as.character(df_melt.pos$variable),"[.][.]"), function(x) list(trimws(x))))
df_melt.pos$family <- sapply(df_melt.pos$levels, function(x) paste(x[5:5], collapse=";"))
df_melt.pos$genus <- sapply(df_melt.pos$levels, function(x) paste(x[1:6], collapse=";"))
df_melt.pos$species <- sapply(df_melt.pos$levels, function(x) paste(x[1:7], collapse=";"))
df_melt.below<-subset(df_melt.pos, family=="f__Enterobacteriaceae")


#only enteros
#GENUS aggregate keeping clusters, patients, and taxa
df.arg.agg.genus <-aggregate(value ~patients + Cluster+ genus, data =  df_melt.below, sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.genus$value <- (df.arg.agg.genus$value>0)+0
colnames(df.arg.agg.genus)<-c("patients","cluster","genus","value")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.genus.reduced <-aggregate(value ~ patients + genus, data =  df.arg.agg.genus, sum)
colnames(df.arg.agg.genus.reduced)<-c("patients","genus","genus_count")

#only enteros
#SPECIES aggregate keeping clusters, patients, and taxa
df.arg.agg.species <-aggregate(value ~patients + Cluster+ species, data =  df_melt.below, sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.species$value <- (df.arg.agg.species$value>0)+0
colnames(df.arg.agg.species)<-c("patients","cluster","species","value")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.species.reduced <-aggregate(value ~ patients + species, data =  df.arg.agg.species, sum)
colnames(df.arg.agg.species.reduced)<-c("patients","species","species_count")


pdf(paste("Unique_ARGs_entero-genera_per_patient_allargs_clustercollapse.pdf",sep=""),height=20, width=10,useDingbats=F)
g <-ggplot(df.arg.agg.genus.reduced, aes(genus,genus_count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(g)
dev.off()

pdf(paste("Unique_ARGs_entero-species_per_patient_allargs_clustercollapse.pdf",sep=""),height=20, width=10,useDingbats=F)
s <-ggplot(df.arg.agg.species.reduced, aes(species,species_count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(s)
dev.off()
################################################################################################
#GENUS aggregate keeping clusters, patients, and taxa
df.arg.agg.genus <-aggregate(value~patients+Cluster+genus,data=df_melt.pos,sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.genus$value <- (df.arg.agg.genus$value>0)+0
colnames(df.arg.agg.genus)<-c("patients","cluster","genus","value")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.genus.reduced <-aggregate(value ~ patients + genus, data =  df.arg.agg.genus, sum)
colnames(df.arg.agg.genus.reduced)<-c("patients","genus","genus_count")

#SPECIES aggregate keeping clusters, patients, and taxa
df.arg.agg.species <-aggregate(value~patients+Cluster+species,data=df_melt.pos,sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.species$value <- (df.arg.agg.species$value>0)+0
colnames(df.arg.agg.species)<-c("patients","cluster","species","value")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.species.reduced <-aggregate(value ~ patients + species, data =  df.arg.agg.species, sum)
colnames(df.arg.agg.species.reduced)<-c("patients","species","species_count")

pdf(paste("Unique_ARGs-genera_per_patient_allargs_clustercollapse.pdf",sep=""),height=20, width=30,useDingbats=F)
g <-ggplot(df.arg.agg.genus.reduced, aes(genus,genus_count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(g)
dev.off()

pdf(paste("Unique_ARGs-species_per_patient_allargs_clustercollapse.pdf",sep=""),height=15, width=40,useDingbats=F)
s <-ggplot(df.arg.agg.species.reduced, aes(species,species_count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(s)
dev.off()
###########################################OXOXOXOXOXOXOXOX####################################
#bring in the bacteria abundances via 
#import max abundance
abund <-read.table("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/all_taxonomic_abundances_basehic_nophage.txt",sep='\t', header = F)
colnames(abund)<-c("Scaffold","Sample","Patient","Abundance","Organism")
organisms = abund$Organism
organisms = gsub("[']", ".",organisms)
organisms = gsub("; ", "..",organisms)
organisms = gsub(";", ".",organisms)
organisms = gsub("[[]", ".",organisms)
organisms = gsub("[]]", ".",organisms)
organisms = gsub("[-]", ".",organisms)
organisms = gsub("[/]", ".",organisms)
organisms = gsub("[ ]", ".",organisms)
abund$Organism <- organisms
#aggregate to total abundance of each organism at family level
abund$Levels <-as.vector(sapply(strsplit(as.character(abund$Organism),"[.][.]"), function(x) list(trimws(x))))
abund$Family <- sapply(abund$Levels, function(x) paste(x[5:5], collapse=";"))
abund.fam.agg <-aggregate(Abundance ~ Sample +  Patient + Family, data = abund, sum)

#get the species
abund$Species <- sapply(abund$Levels, function(x) paste(x[1:7], collapse=";"))
#average all of the contigs based on their sample/patient/species
abund.species.agg <-aggregate(Abundance ~ Sample +  Patient + Species, data = abund, median)
#take the maximal abundance of a particular sample
abund.species.agg.patient <-aggregate(Abundance ~ Patient + Species, data =abund.species.agg, max)

#df.arg.agg.species.reduced$species <- sapply(df.arg.agg.species.reduced$Levels, function(x) paste(x[1:7], collapse=";"))
argcount_speciesabund_merge <- merge(df.arg.agg.species.reduced,abund.species.agg.patient, by.x=c("patients","species"), by.y=c("Patient","Species"),all.x=T)
#"patients"      "species"       "species_count" "Abundance"
argcount_speciesabund_merge$levels = as.vector(sapply(strsplit(as.character(argcount_speciesabund_merge$species),";"), function(x) list(trimws(x))))
argcount_speciesabund_merge$single_species <- sapply(argcount_speciesabund_merge$levels, function(x) paste(x[7:7], collapse=";"))
argcount_speciesabund_merge_speciesonly <- subset(argcount_speciesabund_merge, single_species!="s__.")


fit = lm(species_count~Abundance, data=argcount_speciesabund_merge_speciesonly)
adj.r2 = summary(fit)$adj.r.squared
p = summary(fit)$coefficients[2,4]
pdf(paste("Unique_ARGs-speciesabund_versus_argcount_patient_clustercollapse.pdf",sep=""),height=10, width=10,useDingbats=F)
s <-ggplot(data=argcount_speciesabund_merge_speciesonly, aes(x=Abundance,y=species_count))+
	geom_point(aes(col=patients), position=position_jitter())+
	xlab("Average (Median) Species Abundance")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(s)
dev.off()



###################################################
#top ARGs
###################################################
#another follow-up question
#where are the top-10-ARGs who are they residing in??
top10<-subset(sub5, Top_ARG!="NA")
topgeneral<-read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/top_args_categories.csv", header=T)
top10 <-merge( topgeneral,top10,by="Top_ARG")
top10 <-subset(top10, Confidence!=3)

sub_argcount <- top10
sub_argcount$connectsum <-NULL
sub_argcount$Sample <-NULL
sub_argcount$Confidence <-NULL
sub_argcount$ARG_name<-NULL
sub_argcount$Specific_name<-NULL

df_melt2 <- melt(sub_argcount, id = c("patients", "Cluster", "Top_ARG","Broad_mechanism","Mechanism","Type"))
df_melt.pos <-subset(df_melt2, value>0)
df_melt.pos$levels <-as.vector(sapply(strsplit(as.character(df_melt.pos$variable),"[.][.]"), function(x) list(trimws(x))))
df_melt.pos$family <- sapply(df_melt.pos$levels, function(x) paste(x[1:5], collapse=";"))
df_melt.pos$genus <- sapply(df_melt.pos$levels, function(x) paste(x[1:6], collapse=";"))
df_melt.pos$species <- sapply(df_melt.pos$levels, function(x) paste(x[1:7], collapse=";"))
df_melt.pos$familyname <- sapply(df_melt.pos$levels, function(x) paste(x[5:5], collapse=";"))
df_melt.pos$genusname <- sapply(df_melt.pos$levels, function(x) paste(x[6:6], collapse=";"))
df_melt.pos$speciesname <- sapply(df_melt.pos$levels, function(x) paste(x[6:7], collapse=";"))

#family
df.arg.agg.family <-aggregate(value ~patients+Cluster+Top_ARG+Broad_mechanism+Mechanism+Type+family+familyname, data =  df_melt.pos, sum)
df.arg.agg.family$value <- (df.arg.agg.family$value>0) + 0
colnames(df.arg.agg.family)<-c("Patient","Cluster","Top_ARG","Broad_mechanism","Mechanism","Type","Family","Familyname","Count")
df.arg.agg.family.reduced<-aggregate(Count~Patient+Family + Familyname,data=df.arg.agg.family,sum)
colnames(df.arg.agg.family.reduced)<-c("Patient","Family","Familyname","Count")

#genus
#aggregate keeping clusters, patients, and taxa
df.arg.agg.genus<-aggregate(value~patients+Cluster+Top_ARG+Broad_mechanism+Mechanism+Type+genus+genusname,data=df_melt.pos,sum)
df.arg.agg.genus$value<-(df.arg.agg.genus$value>0) + 0
colnames(df.arg.agg.genus)<-c("Patient","Cluster","Top_ARG","Broad_mechanism","Mechanism","Type","Genus", "Genusname","Count")
df.arg.agg.genus.reduced<-aggregate(Count~Patient+Genus + Genusname,data=df.arg.agg.genus,sum)
colnames(df.arg.agg.genus.reduced)<-c("Patient","Genus","Genusname","Count")

#species
#aggregate keeping clusters, patients, and taxa
df.arg.agg.species <-aggregate(value ~patients+Cluster+Top_ARG+Broad_mechanism+Mechanism+Type+species+speciesname, data =  df_melt.pos, sum)
df.arg.agg.species$value <- (df.arg.agg.species$value>0) + 0
colnames(df.arg.agg.species)<-c("Patient","Cluster","Top_ARG","Broad_mechanism","Mechanism","Type","Species", "Speciesname","Count")
df.arg.agg.species.reduced<-aggregate(Count~Patient+Species+Speciesname,data=df.arg.agg.species,sum)
colnames(df.arg.agg.species.reduced)<-c("Patient","Species","Speciesname","Count")



#################################################################

pdf(paste("Unique_TopARGs_family_per_patient_clustercollapse.pdf",sep=""),height=10, width=10,useDingbats=F)
f <-ggplot(df.arg.agg.family.reduced, aes(Family,Count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=Patient), position=position_jitter(w = 0.2, h = 0))+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(f)
dev.off()

pdf(paste("Unique_TopARGs_genera_per_patient_clustercollapse.pdf",sep=""),height=10, width=10,useDingbats=F)
g <-ggplot(df.arg.agg.genus.reduced, aes(Genus,Count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=Patient), position=position_jitter(w = 0.2, h = 0))+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(g)
dev.off()

pdf(paste("Unique_TopARGs_species_per_patient_clustercollapse.pdf",sep=""),height=15, width=20,useDingbats=F)
s <-ggplot(df.arg.agg.species.reduced, aes(Species,Count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=Patient), position=position_jitter(w = 0.2, h = 0))+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(s)
dev.off()


##############top arg categories vs family
pdf(paste("Unique_TopARGs_categories_family_bypatient_clustercollapse.pdf",sep=""),height=15, width=20,useDingbats=F)
	f <-ggplot(df.arg.agg.family, aes(x = Family, y = Count, fill = Broad_mechanism))+
	geom_bar(stat="identity")+
	facet_grid(.~Patient)+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(f)
dev.off()

pdf(paste("Unique_TopARGs_categories_family_clustercollapse.pdf",sep=""),height=15, width=8,useDingbats=F)
	f <-ggplot(df.arg.agg.family, aes(x = Family, y = Count, fill = Broad_mechanism))+
	geom_bar(stat="identity")+
	xlab("Taxa")+
	ylab("ARG Count--Patients are unique ARGs")+
	scale_color_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(f)
dev.off()

################################################################
#this will collapse on the patients and keep the clusters unique
df.arg.agg.family.patient.reduced<-aggregate(Count~Family+Familyname+Cluster+Broad_mechanism,data=df.arg.agg.family,sum)
df.arg.agg.family.patient.reduced$Count <- (df.arg.agg.family.patient.reduced$Count>0) + 0
colnames(df.arg.agg.family.patient.reduced)<-c("Family","Familyname","Cluster","Broad_mechanism","Count")
pdf(paste("Unique_TopARGs_categories_family_clusterpatientcollapse.pdf",sep=""),height=20, width=15,useDingbats=F)
	f <-ggplot(df.arg.agg.family.patient.reduced, aes(x = Family, y = Count, fill = Broad_mechanism))+
	geom_bar(stat="identity")+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(f)
dev.off()

#this will collapse on the patients and keep the clusters unique
df.arg.agg.genus.patient.reduced<-aggregate(Count~Genus+Genusname+Cluster+Broad_mechanism,data=df.arg.agg.genus,sum)
df.arg.agg.genus.patient.reduced$Count <- (df.arg.agg.genus.patient.reduced$Count>0) + 0
colnames(df.arg.agg.genus.patient.reduced)<-c("Genus","Genusname","Cluster","Broad_mechanism","Count")
pdf(paste("Unique_TopARGs_categories_genus_clusterpatientcollapse.pdf",sep=""),height=20, width=15,useDingbats=F)
	f <-ggplot(df.arg.agg.genus.patient.reduced, aes(x = Genus, y = Count, fill = Broad_mechanism))+
	geom_bar(stat="identity")+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(f)
dev.off()
######################################################################################
df.arg.agg.species.patient.reduced<-aggregate(Count~Species+Speciesname+Cluster+Broad_mechanism,data=df.arg.agg.species,sum)
df.arg.agg.species.patient.reduced$Count <- (df.arg.agg.species.patient.reduced$Count>0) + 0
colnames(df.arg.agg.species.patient.reduced)<-c("Species","Speciesname","Cluster","Broad_mechanism","Count")
pdf(paste("Unique_TopARGs_categories_species_clusterpatientcollapse.pdf",sep=""),height=20, width=20,useDingbats=F)
	f <-ggplot(df.arg.agg.species.patient.reduced, aes(x = Species, y = Count, fill = Broad_mechanism))+
	geom_bar(stat="identity")+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(f)
dev.off()

###############################so can we ask...
#plot species abundances versus ARG abundances
pdf("TopARGS_ARG_vs_Species_abundances.pdf")
species.argrpkm <- merge(df.arg.agg.species, df5.agg, by=c("Patient","Cluster"))
argrpkm.speciesabundance<-merge(species.argrpkm, abund.species.agg.patient, by=c("Patient","Species"))
	f <-ggplot(argrpkm.speciesabundance, aes(x = Abundance, y = RPKM))+
	geom_point(aes(col=Patient))+
	xlab("Taxa Abundance")+
	ylab("ARG Abundance (RPKM)")+
	scale_color_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
	print(f)
dev.off()





