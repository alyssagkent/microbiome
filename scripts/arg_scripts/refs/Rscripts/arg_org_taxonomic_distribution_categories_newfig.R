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
#df_melt.below<-subset(df_melt.pos, family!="f__")


#aggregate keeping clusters, patients, and taxa
df.arg.agg <-aggregate(value ~ patients +  Cluster+variable + family + familyname, data=df_melt.pos, FUN=sum)
df.arg.agg$value <- (df.arg.agg$value>0)+0
colnames(df.arg.agg)<-c("patients","cluster","variable","family","familyname","value")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.fam <-aggregate(value ~ cluster + patients +family+familyname, data =  df.arg.agg, sum)
colnames(df.arg.agg.fam)<-c("cluster","patients","family","familyname","value")
df.arg.agg.fam$value <- (df.arg.agg.fam$value>0)+0
df.arg.agg.fam.arginfo <- merge(arginfo, df.arg.agg.fam, by.x="Cluster", by.y="cluster", all.y=T)
#se
df.arg.agg.fam.arginfo.reduced <-aggregate(value ~ Sub_mechanism +  patients +family+familyname, data =  df.arg.agg.fam.arginfo, sum)
df.arg.agg.fam.arginfo.reduced.agg<-aggregate(value ~patients +family+familyname, data =  df.arg.agg.fam.arginfo, sum)


pdf(paste("Unique_ARGs_family_patient_categories_allargs_clustercollapse.pdf",sep=""),height=30, width=10,useDingbats=F)
plot1 <-ggplot(data=df.arg.agg.fam.arginfo.reduced, aes(factor(familyname),value))+
	geom_boxplot(outlier.shape=NA)+
	facet_grid(Sub_mechanism~.,scales="free_y")+ 
	geom_point(aes(col=patients), position=position_jitter(height = 0))+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),strip.text.y = element_text(angle = 0))
print(plot1)
dev.off()


pdf(paste("Unique_ARGs_family_patient_categories_allargs_clustercollapse_barplot.pdf",sep=""),height=60, width=10,useDingbats=F)
plot2 <-ggplot(data=df.arg.agg.fam.arginfo.reduced, aes(factor(familyname),value))+
        geom_bar(aes(fill=patients),stat="identity")+
        facet_grid(Sub_mechanism~.,scales="free_y") +
        xlab("Taxa")+
        ylab("ARG Count")+
        scale_fill_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot2)
dev.off()

pdf(paste("Unique_ARGs_family_patient_allargs_clustercollapse.pdf",sep=""),height=10, width=10,useDingbats=F)
plot3 <-ggplot(data=df.arg.agg.fam.arginfo.reduced.agg, aes(factor(familyname),value))+
        geom_boxplot(outlier.shape=NA)+
        geom_point(aes(col=patients), position=position_jitter(height = 0))+
        xlab("Taxa")+
        ylab("ARG Count")+
        scale_color_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),strip.text.y = element_text(angle = 0))
print(plot3)
dev.off()

################################################################################################
#GENUS aggregate keeping clusters, patients, and taxa
df_melt.pos$genus <- sapply(df_melt.pos$levels, function(x) paste(x[1:6], collapse=";"))
df_melt.pos$genusname <- sapply(df_melt.pos$levels, function(x) paste(x[6:6], collapse=";"))
df_melt.pos$species <- sapply(df_melt.pos$levels, function(x) paste(x[1:7], collapse=";"))
df_melt.pos$speciesname <- sapply(df_melt.pos$levels, function(x) paste(x[7:7], collapse=";"))

df.arg.agg.genus <-aggregate(value~patients+Cluster+family+genus+genusname,data=df_melt.pos,sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.genus$value <- (df.arg.agg.genus$value>0)+0
colnames(df.arg.agg.genus)<-c("Patient","Cluster","Family","Genus","Genusname","Count")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.genus.reduced <-aggregate(Count ~ Patient + Family + Genus+Genusname, data =  df.arg.agg.genus, sum)
colnames(df.arg.agg.genus.reduced)<-c("Patient","Family","Genus","Genusname","Genus_count")

#SPECIES aggregate keeping clusters, patients, and taxa
df.arg.agg.species <-aggregate(value~patients+Cluster+family+species+speciesname,data=df_melt.pos,sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.species$value <- (df.arg.agg.species$value>0)+0
colnames(df.arg.agg.species)<-c("Patient","Cluster","Family","Species","Speciesname","Count")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.species.reduced <-aggregate(Count ~ Patient + Family+Species+Speciesname, data =  df.arg.agg.species, sum)
colnames(df.arg.agg.species.reduced)<-c("Patient","Family","Species","Speciesname","Species_count")

pdf(paste("Unique_ARGs-genera_per_patient_allargs_clustercollapse.pdf",sep=""),height=20, width=30,useDingbats=F)
plot4 <-ggplot(df.arg.agg.genus.reduced, aes(Genus,Genus_count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=Patient), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot4)
dev.off()

pdf(paste("Unique_ARGs-species_per_patient_allargs_clustercollapse.pdf",sep=""),height=15, width=40,useDingbats=F)
plot5 <-ggplot(df.arg.agg.species.reduced, aes(Species,Species_count))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=Patient), position=position_jitter())+
	xlab("Taxa")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot5)
dev.off()

###############################trim to enteros
df.arg.agg.genus.reduced.entero <- subset(df.arg.agg.genus.reduced, Family == "f__Enterobacteriaceae")
pdf(paste("Unique_ARGs-genera_per_patient_allargs_clustercollapse_enteros.pdf",sep=""),height=20, width=30,useDingbats=F)
plot4_entero <-ggplot(df.arg.agg.genus.reduced.entero, aes(Genus,Genus_count))+
        geom_boxplot(outlier.shape=NA)+
        geom_point(aes(col=Patient), position=position_jitter())+
        xlab("Taxa")+
        ylab("ARG Count")+
        scale_color_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot4_entero)
dev.off()

df.arg.agg.species.reduced.entero <- subset(df.arg.agg.species.reduced, Family == "f__Enterobacteriaceae")
pdf(paste("Unique_ARGs-species_per_patient_allargs_clustercollapse_enteros.pdf",sep=""),height=15, width=40,useDingbats=F)
plot5_entero <-ggplot(df.arg.agg.species.reduced.entero, aes(Species,Species_count))+
        geom_boxplot(outlier.shape=NA)+
        geom_point(aes(col=Patient), position=position_jitter())+
        xlab("Taxa")+
        ylab("ARG Count")+
        scale_color_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot5_entero)
dev.off()
################################################################
#ok now merge it with the arg abundances for each person??
sub_argcount <- sub5
sub_argcount$connectsum <- NULL
sub_argcount$Top_ARG<-NULL
sub_argcount$ARG_name<-as.character(sub_argcount$ARG_name)
df_melt2 <- melt(sub_argcount, id = c("Sample","patients", "Cluster", "ARG_name"))
df_melt.pos <-subset(df_melt2, value>0)
df_melt.pos$levels <-as.vector(sapply(strsplit(as.character(df_melt.pos$variable),"[.][.]"), function(x) list(trimws(x))))
df_melt.pos$family <- sapply(df_melt.pos$levels, function(x) paste(x[5:5], collapse=";"))
df_melt.pos$familyname <- sapply(df_melt.pos$levels, function(x) paste(x[1:5], collapse=";"))
df_melt.pos$species <- sapply(df_melt.pos$levels, function(x) paste(x[1:7], collapse=";"))
df_melt.pos$speciesname <- sapply(df_melt.pos$levels, function(x) paste(x[7:7], collapse=";"))

df.arg.agg.species <-aggregate(value~Sample + patients+Cluster+family+species+speciesname,data=df_melt.pos,sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.species$value <- (df.arg.agg.species$value>0)+0
colnames(df.arg.agg.species)<-c("Sample","Patient","Cluster","Family","Species","Speciesname","Count")
df.arg.agg.species.argrpkm <- merge(df.arg.agg.species, df5.melt, by=c("Sample","Patient","Cluster"))
df.arg.agg.species.reduced <-aggregate(RPKM ~ Patient + Family+Species+Speciesname, data =  df.arg.agg.species.argrpkm, sum)
colnames(df.arg.agg.species.reduced)<-c("Patient","Family","Species","Speciesname","RPKM")
df.arg.agg.species.reduced.entero <- subset(df.arg.agg.species.reduced, Family == "f__Enterobacteriaceae")

pdf("Unique_ARGs-species_per_timepoint_allargs_clustercollapse_enteros_argrpkmssummed.pdf",height=20, width=30,useDingbats=F)
ggplot(df.arg.agg.species.reduced.entero, aes(Species, RPKM))+
	geom_boxplot(outlier.shape=NA)+
	geom_point(aes(col=Patient),position = position_jitter())+
	xlab("Taxa")+
	ylab("ARG Summed Abundance")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

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
abund$Levels <- as.vector(sapply(strsplit(as.character(abund$Organism),"[.][.]"), function(x) list(trimws(x))))

#get the species
abund$Species <- sapply(abund$Levels, function(x) paste(x[1:7], collapse=";"))
#average all of the contigs based on their sample/patient/species
abund.species.agg <-aggregate(Abundance ~ Sample +  Patient + Species, data = abund, median)
#take the maximal abundance of a particular sample
abund.species.agg.patient <-aggregate(Abundance ~ Patient + Species, data =abund.species.agg, max)

#aggregate to total abundance of each organism at family level
abund.species.agg.patient$Levels <-as.vector(sapply(strsplit(as.character(abund.species.agg.patient$Species),";"), function(x) list(trimws(x))))
abund.species.agg.patient$Family <- sapply(abund.species.agg.patient$Levels, function(x) paste(x[1:5], collapse=";"))
abund.spp.agg.pat.fam <-aggregate(Abundance ~ Patient + Family, data = abund.species.agg.patient, sum)

#just plot the abundance of each of the families
df.argfam.abund<-merge(df.arg.agg.fam.arginfo.reduced, abund.spp.agg.pat.fam, by.x=c("familyname", "patients"), by.y=c("Family","Patient"), all.x=T)
df.argfam.abund$Sub_mechanism <-NULL
df.argfam.abund$Mechanism<-NULL
df.argfam.abund$value <-NULL
df.famwitharg.abund <- unique(df.argfam.abund)
pdf("Family_abundances_sum_max_median_barplot.pdf",height=15, width=10,useDingbats=F)
plot6 <- ggplot(data=df.famwitharg.abund, aes(familyname, Abundance))+
        geom_bar(aes(fill=patients),stat="identity")+
        xlab("Taxa")+
        ylab("Taxa Abundance--Sum(Max(Median(Species Contigs)))")+
        scale_fill_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot6)
dev.off()

pdf("Family_abundances_sum_max_median_boxplot.pdf",height=15, width=10,useDingbats=F)
plot6_2 <- ggplot(data=df.famwitharg.abund, aes(familyname, Abundance))+
        geom_boxplot(outlier.shape=NA)+
        geom_point(aes(col=patients))+
        xlab("Taxa")+
        ylab("Taxa Abundance--Sum(Max(Median(Species Contigs)))")+
        scale_colour_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot6_2)
dev.off()

#df.arg.agg.species.reduced$species <- sapply(df.arg.agg.species.reduced$Levels, function(x) paste(x[1:7], collapse=";"))
argcount_speciesabund_merge <- merge(df.arg.agg.species.reduced,abund.species.agg.patient, by=c("Patient","Species"),all.x=T)
argcount_speciesabund_merge$Levels = as.vector(sapply(strsplit(as.character(argcount_speciesabund_merge$Species),";"), function(x) list(trimws(x))))
argcount_speciesabund_merge$Single_species <- sapply(argcount_speciesabund_merge$Levels, function(x) paste(x[7:7], collapse=";"))
argcount_speciesabund_merge_speciesonly <- subset(argcount_speciesabund_merge, Single_species!="s__.")

fit = lm(Species_count~Abundance, data=argcount_speciesabund_merge_speciesonly)
adj.r2 = summary(fit)$adj.r.squared
p = summary(fit)$coefficients[2,4]
pdf(paste("Unique_ARGs-speciesabund_versus_argcount_patient_clustercollapse.pdf",sep=""),height=10, width=10,useDingbats=F)
plot7 <-ggplot(data=argcount_speciesabund_merge_speciesonly, aes(x=Abundance,y=Species_count))+
	geom_point(aes(col=Patient), position=position_jitter())+
	xlab("Max of Timepoints(Median(Species Contig RPKM))")+
	ylab("ARG Count")+
	scale_color_manual(label=pats, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot7)
dev.off()



####################################################################################################
#lasso plot with ALL AR genes
pdf("Allargs_vs_Species_abundances_lasso.pdf",height=7,width=9.5,useDingbats=F)
species.argrpkm <- merge(df.arg.agg.species, df5.agg, by=c("Patient","Cluster"))
onlyspecies.argrpkm <- subset(species.argrpkm, Speciesname!="s__.")
argrpkm.speciesabundance<-merge(onlyspecies.argrpkm, abund.species.agg.patient, by=c("Patient","Species"))
argrpkm.speciesabund.lasso<-aggregate(Abundance~Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df.to.plot <- subset(argrpkm.speciesabund.lasso)
plot8 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        xlim(0,430)+
        ylim(0,430)+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+ 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot8)
fit = lm(df.to.plot$RPKM~df.to.plot$Abundance)
summary(fit)
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="pearson")
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="spearman")
dev.off()

pdf("Allargs_vs_Species_abundances_lasso_wo_outlier.pdf",height=7,width=9.5,useDingbats=F)
plot9 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
	xlim(0,430)+
	ylim(0,170)+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
	scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot9)
dev.off()

pdf(paste("Allargs_vs_Species_abundances_lasso_separated.pdf",sep=""),height=60,width=9.5,useDingbats=F)
argrpkm.speciesabund.lasso<-aggregate(Abundance~Patient+Cluster+RPKM+Sub_mechanism,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Patient+Cluster+RPKM+Sub_mechanism,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df.to.plot <- subset(argrpkm.speciesabund.lasso)
plot10 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
	facet_grid(Sub_mechanism~.)+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot10)
dev.off()

pdf(paste("Allargs_vs_Species_abundances_lasso_separated_nooutlier.pdf",sep=""),height=60,width=9.5,useDingbats=F)
plot11 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        ylim(0,170)+
        facet_grid(Sub_mechanism~.)+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot11)
dev.off()


#########################################################
#lassoplots connection was ever present
#for each cluster-patient-speciesgroup combo from above...can you get the abundances at each timepoint for that comb
onlyspecies.arg <- subset(df.arg.agg.species, Speciesname!="s__.")
#bring in somehow all the timepoints before you aggregate
#ugh this is missing the zeroes for timepoints
#first merge abundance data with a file that has all timepoints vs all species
samples = unique(abund.species.agg$Sample)
species = unique(abund.species.agg$Species)
allcombinations = expand.grid(samples,species)
colnames(allcombinations)<-c("Sample","Species")
abund.species.agg.allcombos <- merge(allcombinations, abund.species.agg, by=c("Sample","Species"),all.x=T)
abund.species.agg.allcombos$Patient<-as.vector(sapply(strsplit(as.character(abund.species.agg.allcombos$Sample),"-"), function(x) x[[1]]))
abund.species.agg.allcombos$Abundance[is.na(abund.species.agg.allcombos$Abundance)]<-0
arg.speciesabundance<-merge(onlyspecies.arg, abund.species.agg.allcombos, by=c("Patient","Species"),all.x=T)

#bring in the ARG abundances based on the sample
argrpkm.speciesabundance <- merge(arg.speciesabundance, df5.melt, by=c("Sample","Patient","Cluster"))
argrpkm.speciesabund.lasso<-aggregate(Abundance~Sample +Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Sample+Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df.to.plot <- subset(argrpkm.speciesabund.lasso)
pdf("Allargs_vs_Species_abundances_lasso_samples-presentever.pdf",height=7,width=9.5,useDingbats=F)
plot12 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot12)
fit = lm(df.to.plot$RPKM~df.to.plot$Abundance)
summary(fit)
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="pearson")
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="spearman")
dev.off()

pdf("Allargs_vs_Species_abundances_lasso_samples-presentever_nooutlier.pdf",height=7,width=9.5,useDingbats=F)
plot13 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
	ylim(0,170)+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot13)
dev.off()

pdf(paste("Allargs_vs_Species_abundances_lasso_samples-presentever_separated.pdf",sep=""),height=60,width=9.5,useDingbats=F)
argrpkm.speciesabund.lasso<-aggregate(Abundance~Sample +Patient+Cluster+RPKM+Sub_mechanism,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Sample+Patient+Cluster+RPKM+Sub_mechanism,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df.to.plot <- subset(argrpkm.speciesabund.lasso, Count<50) 
plot14 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        facet_grid(Sub_mechanism~.)+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot14)
dev.off()

pdf(paste("Allargs_vs_Species_abundances_lasso_samples-presentever_separated_nooutlier.pdf",sep=""),height=60,width=9.5,useDingbats=F)
plot15 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        ylim(0,170)+
        facet_grid(Sub_mechanism~.)+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot15)
dev.off()


#only with the connections that are present in the sample----spotty connections undermine this
sub7 <- sub5
sub7$Top_ARG<-NULL
sub7$ARG_name<-as.character(sub7$ARG_name)
sub_argcount <- sub7
sub_argcount$connectsum <- NULL
df_melt2 <- melt(sub_argcount, id = c("Sample","patients", "Cluster", "ARG_name"))
df_melt.pos <-subset(df_melt2, value>0)
df_melt.pos$levels <-as.vector(sapply(strsplit(as.character(df_melt.pos$variable),"[.][.]"), function(x) list(trimws(x))))
df_melt.pos$species <- sapply(df_melt.pos$levels, function(x) paste(x[1:7], collapse=";"))
df_melt.pos$speciesname <- sapply(df_melt.pos$levels, function(x) paste(x[7:7], collapse=";"))
#SPECIES aggregate keeping clusters, patients, and taxa
df.arg.agg.species <-aggregate(value~Sample+patients+Cluster+species+speciesname,data=df_melt.pos,sum)
#make binary so patients aren't contributing more because of time-points
df.arg.agg.species$value <- (df.arg.agg.species$value>0)+0
colnames(df.arg.agg.species)<-c("Sample","Patient","Cluster","Species","Speciesname","Count")

pdf("Allargs_vs_Species_abundances_lasso_samples-present.pdf",height=7,width=9.5,useDingbats=F)
species.argrpkm <- merge(df.arg.agg.species, df5.melt, by=c("Sample","Patient","Cluster"))
onlyspecies.argrpkm <- subset(species.argrpkm, Speciesname!="s__.")
argrpkm.speciesabundance<-merge(onlyspecies.argrpkm, abund.species.agg, by=c("Sample","Patient","Species"))
argrpkm.speciesabund.lasso<-aggregate(Abundance~Sample+Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Sample+Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df.to.plot <- subset(argrpkm.speciesabund.lasso, Count<50)
plot16 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot16)
fit = lm(df.to.plot$RPKM~df.to.plot$Abundance)
summary(fit)
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="pearson")
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="spearman")
dev.off()


pdf(paste("Allargs_vs_Species_abundances_lasso_samples-present_separated.pdf",sep=""),height=60,width=9.5,useDingbats=F)
argrpkm.speciesabund.lasso<-aggregate(Abundance~Sample+Patient+Cluster+RPKM+Sub_mechanism,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Sample+Patient+Cluster+RPKM+Sub_mechanism,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df.to.plot <- subset(argrpkm.speciesabund.lasso)
plot17 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient, size = Count))+
        facet_grid(Sub_mechanism~.)+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot17)
dev.off()


#make a table of the correlations
s1 = subset(argrpkm.speciesabund.lasso, Sub_mechanism=="Acetyltransferase")
print(unique(s1$Sub_mechanism))
cor.test(s1$Abundance, s1$RPKM)
s1 = subset(argrpkm.speciesabund.lasso, Sub_mechanism=="Beta-Lactamase")
print(unique(s1$Sub_mechanism))
cor.test(s1$Abundance, s1$RPKM)
s1 = subset(argrpkm.speciesabund.lasso, Sub_mechanism=="Efflux")
print(unique(s1$Sub_mechanism))
cor.test(s1$Abundance, s1$RPKM)
s1 = subset(argrpkm.speciesabund.lasso, Sub_mechanism=="Gene Modulating Resistance")
print(unique(s1$Sub_mechanism))
cor.test(s1$Abundance, s1$RPKM)
s1 = subset(argrpkm.speciesabund.lasso, Sub_mechanism=="Quinolone Resistance")
print(unique(s1$Sub_mechanism))
cor.test(s1$Abundance, s1$RPKM)
s1 = subset(argrpkm.speciesabund.lasso, Sub_mechanism=="Tetracycline Resistance")
print(unique(s1$Sub_mechanism))
cor.test(s1$Abundance, s1$RPKM,method="pearson")
cor.test(s1$Abundance, s1$RPKM,method="spearman")


df <- data.frame(matrix(ncol = 15, nrow = 14))
i = 1
for (mech_class in unique(argrpkm.speciesabund.lasso$Sub_mechanism)){
s1 = subset(argrpkm.speciesabund.lasso, Sub_mechanism==mech_class)
print(mech_class)
a = cor.test(s1$Abundance, s1$RPKM,method="pearson")
b = cor.test(s1$Abundance, s1$RPKM,method="spearman")
c = lm(RPKM~Abundance, data=s1)
d = summary(c)
val = c(mech_class, a$estimate, a$p.value, a$t, a$parameter, b$estimate, b$p.value, d$coefficients[2,], d$fstatistic, d$r.squared, d$adj.r.squared)
df[i,] = val
i = i + 1
}

colnames(df)<-c("Mechanism","Pearson_r","Pearson_p-value","DF","Spearman_rho","Spearman_p-value","LinReg_Abundance_estimate","LinReg_Std.Error","LinReg_t-value","LinReg_p-value","LinReg_fstatistic","LinReg_fstatistic_df1","LinReg_fstatistic_df2","LinReg_R-sq","LinReg_Adj.R-sq")

write.table(df, "Allargs_vs_Species_abundances_lasso_samples-present_separated_statistics.txt", sep="\t",row.names=F)



#############################################################################sample family abundances

abund.species.agg$Levels <-as.vector(sapply(strsplit(as.character(abund.species.agg$Species),";"), function(x) list(trimws(x))))
abund.species.agg$Family <- sapply(abund.species.agg$Levels, function(x) paste(x[1:5], collapse=";"))
abund.spp.agg.fam <-aggregate(Abundance ~ Sample +Patient + Family, data = abund.species.agg, sum)

#just plot the abundance of each of the families
#aggregate keeping clusters, patients, and taxa
df_melt.pos$family <- sapply(df_melt.pos$levels, function(x) paste(x[1:5], collapse=";"))
df_melt.pos$familyname <- sapply(df_melt.pos$levels, function(x) paste(x[5:5], collapse=";"))
df.arg.agg <-aggregate(value ~ Sample+patients +  Cluster+variable + family + familyname, data=df_melt.pos, FUN=sum)
df.arg.agg$Count <- (df.arg.agg$value>0)+0
colnames(df.arg.agg)<-c("Sample","Patient","Cluster","variable","Family","Familyname","Count")
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.fam <-aggregate(Count ~ Cluster+Sample + Patient +Family+Familyname, data =  df.arg.agg, sum)
colnames(df.arg.agg.fam)<-c("Cluster","Sample","Patient","Family","Familyname","Count")
df.arg.agg.fam$Count <- (df.arg.agg.fam$Count>0)+0
df.arg.agg.fam.arginfo <- merge(arginfo, df.arg.agg.fam, by="Cluster", all.y=T)
#se
df.arg.agg.fam.arginfo.reduced <-aggregate(Count ~ Sub_mechanism + Mechanism + Sample+ Patient +Family+Familyname, data =  df.arg.agg.fam.arginfo, sum)
df.arg.agg.fam.arginfo.reduced.agg<-aggregate(Count ~Sample+Patient +Family+Familyname, data =  df.arg.agg.fam.arginfo, sum)
df.argfam.abund<-merge(df.arg.agg.fam.arginfo.reduced, abund.spp.agg.fam, by=c("Family","Sample","Patient"), all.x=T)
df.argfam.abund$Sub_mechanism <-NULL
df.argfam.abund$Mechanism<-NULL
df.argfam.abund$Count <-NULL
df.famwitharg.abund <- unique(df.argfam.abund)
pdf("Family_abundances_samples_boxplot.pdf",height=15, width=10,useDingbats=F)
plot17_2 <- ggplot(data=df.famwitharg.abund, aes(Family, Abundance))+
        geom_boxplot(outlier.shape=NA)+
        geom_point(aes(col=Patient))+
        xlab("Taxa")+
        ylab("Taxa Abundance--Sum((Median(Species Contigs))")+
        scale_colour_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot17_2)
dev.off()


#############################################################SAMPLES
#the problem is that you aren't keeping the samples separate and you're using the weird max abundance thing
#aggregate keeping  patients, and taxa and summing over the args
df.arg.agg.species.reduced <-aggregate(Count ~Sample+ Patient + Species+Speciesname, data =  df.arg.agg.species, sum)
colnames(df.arg.agg.species.reduced)<-c("Sample","Patient","Species","Speciesname","Species_count")
argcount_speciesabund_merge <- merge(df.arg.agg.species.reduced,abund.species.agg, by=c("Sample","Patient","Species"),all.x=T)
argcount_speciesabund_merge$Levels = as.vector(sapply(strsplit(as.character(argcount_speciesabund_merge$Species),";"), function(x) list(trimws(x))))
argcount_speciesabund_merge_speciesonly <- subset(argcount_speciesabund_merge, Speciesname!="s__.")
fit = lm(Species_count~Abundance, data=argcount_speciesabund_merge_speciesonly)
adj.r2 = summary(fit)$adj.r.squared
p = summary(fit)$coefficients[2,4]
pdf(paste("Unique_ARGs-speciesabund_versus_argcount_samples-present_clustercollapse.pdf",sep=""),height=10, width=10,useDingbats=F)
plot18 <-ggplot(data=argcount_speciesabund_merge_speciesonly, aes(x=Abundance,y=Species_count))+
        geom_point(aes(col=Patient), position=position_jitter())+
        xlab("Median(Species Contig RPKM)")+
        ylab("ARG Count")+
        scale_color_manual(label=pats, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot18)
dev.off()


#######################################################################################
#######################################################################################

###################################################
#top ARGs
###################################################
#another follow-up question
#where are the top-10-ARGs who are they residing in??
top10<-subset(sub5, Top_ARG!="NA")
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

#genus
#aggregate keeping clusters, patients, and taxa
df.arg.agg.genus<-aggregate(value~patients+Cluster+Mergedname+Category+familyname+genus+genusname,data=df_melt.pos,sum)
df.arg.agg.genus$value<-(df.arg.agg.genus$value>0) + 0
colnames(df.arg.agg.genus)<-c("Patient","Cluster","Mergedname","Category","Family","Genus","Genusname","Count")
df.arg.agg.genus.reduced<-aggregate(Count~Patient+Family +Genus + Genusname,data=df.arg.agg.genus,sum)
colnames(df.arg.agg.genus.reduced)<-c("Patient","Family","Genus","Genusname","Count")

#species
#aggregate keeping clusters, patients, and taxa
df.arg.agg.species <-aggregate(value ~patients+Cluster+Mergedname+Category+familyname +species+speciesname, data =  df_melt.pos, sum)
df.arg.agg.species$value <- (df.arg.agg.species$value>0) + 0
colnames(df.arg.agg.species)<-c("Patient","Cluster","Mergedname","Category","Family","Species", "Speciesname","Count")
df.arg.agg.species.reduced<-aggregate(Count~Patient+Family+Species+Speciesname,data=df.arg.agg.species,sum)
colnames(df.arg.agg.species.reduced)<-c("Patient","Family","Species","Speciesname","Count")

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

###################################ok do something similar but for the species and the genus---only enteros?

df.arg.agg.species.patient.reduced<-aggregate(Count~Family+Species+Speciesname+Cluster+Mergedname+Category,data=df.arg.agg.species,sum)
df.arg.agg.species.patient.reduced$Count <- (df.arg.agg.species.patient.reduced$Count>0) + 0
colnames(df.arg.agg.species.patient.reduced)<-c("Family","Species","Speciesname","Cluster","Mergedname","Category","Count")
pdf(paste("Unique_TopARGs_categories_species_clusterpatientcollapse.pdf",sep=""),height=20, width=15,useDingbats=F)
plot19_species <-ggplot(df.arg.agg.species.patient.reduced, aes(x = Species, y = Count, fill = Category))+
        geom_bar(stat="identity")+
        xlab("Taxa")+
        ylab("ARG Count")+
        #scale_fill_manual(label=patients, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot19_species)
dev.off()

#subset to only include enteros
df.arg.agg.species.patient.reduced.enteros <- subset(df.arg.agg.species.patient.reduced, Family =="f__Enterobacteriaceae")
pdf(paste("Unique_TopARGs_categories_species_clusterpatientcollapse_enteros.pdf",sep=""),height=20, width=15,useDingbats=F)
plot19_species_enteros <-ggplot(df.arg.agg.species.patient.reduced.enteros, aes(x = Species, y = Count, fill = Category))+
        geom_bar(stat="identity")+
        xlab("Taxa")+
        ylab("ARG Count")+
        #scale_fill_manual(label=patients, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot19_species_enteros)
dev.off()

df.arg.agg.genus.patient.reduced<-aggregate(Count~Family+Genus+Genusname+Cluster+Mergedname+Category,data=df.arg.agg.genus,sum)
df.arg.agg.genus.patient.reduced$Count <- (df.arg.agg.genus.patient.reduced$Count>0) + 0
colnames(df.arg.agg.genus.patient.reduced)<-c("Family","Genus","Genusname","Cluster","Mergedname","Category","Count")
pdf(paste("Unique_TopARGs_categories_genus_clusterpatientcollapse.pdf",sep=""),height=20, width=15,useDingbats=F)
plot19_genus <-ggplot(df.arg.agg.genus.patient.reduced, aes(x = Genus, y = Count, fill = Category))+
        geom_bar(stat="identity")+
        xlab("Taxa")+
        ylab("ARG Count")+
        #scale_fill_manual(label=patients, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot19_genus)
dev.off()

#subset to only include enteros
df.arg.agg.genus.patient.reduced.enteros <- subset(df.arg.agg.genus.patient.reduced, Family =="f__Enterobacteriaceae")
pdf(paste("Unique_TopARGs_categories_genus_clusterpatientcollapse_enteros.pdf",sep=""),height=20, width=15,useDingbats=F)
plot19_genus_enteros <-ggplot(df.arg.agg.genus.patient.reduced.enteros, aes(x = Genus, y = Count, fill = Category))+
        geom_bar(stat="identity")+
        xlab("Taxa")+
        ylab("ARG Count")+
        #scale_fill_manual(label=patients, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot19_genus_enteros)
dev.off()
        



dcolnames(f.arg.agg.genus)<-c("Patient","Cluster","Top_ARG","Broad_mechanism","Mechanism","Type","Genus", "Genusname","Count")

###############################so can we ask...
#plot species abundances versus ARG abundances
pdf("TopARGS_vs_Species_abundances_patients.pdf",useDingbats=F)
species.argrpkm <- merge(df.arg.agg.species, df5.agg, by=c("Patient","Cluster"))
argrpkm.speciesabundance<-merge(species.argrpkm, abund.species.agg.patient, by=c("Patient","Species"))
plot20 <-ggplot(argrpkm.speciesabundance, aes(x = Abundance, y = RPKM))+
	geom_point(aes(col=Patient))+
	xlab("Taxa Abundance")+
	ylab("ARG Abundance (RPKM)")+
	scale_color_manual(label=patients, values=patient.col)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot20)
dev.off()


#next goal is to take these top args and add all of the taxa that represent them together...so that you only have one line for each arg...witha summed species abundance
#is it a straight line? then the species are representative of total populations
#is it not a straight line?


pdf("TopARGS_vs_Species_abundances_lasso.pdf",height=7, width=9.5,useDingbats=F)
species.argrpkm <- merge(df.arg.agg.species, df5.agg, by=c("Patient","Cluster"))
onlyspecies.argrpkm <- subset(species.argrpkm, Speciesname!="s__.")
argrpkm.speciesabundance<-merge(onlyspecies.argrpkm, abund.species.agg.patient, by=c("Patient","Species"))
argrpkm.speciesabund.lasso<-aggregate(Abundance~Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df <- merge(argrpkm.speciesabund.lasso, toparg_clusters, by="Cluster", all.x=T)
df <-merge(known_cat,df,by="Mergedname")
df.to.plot<- df
plot21 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_text(aes(label=Trim_name, cex=2),hjust=-.1, vjust=0)+
        geom_point(aes(col=Patient, size = Count))+
        xlab("Species Abundance (Max of Patient Timepoints(Median(Species Contig RPKMs)))")+
        ylab("ARG Abundance Max of Patient Timepoints(RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot21)
fit <-lm(argrpkm.speciesabund.lasso$RPKM~argrpkm.speciesabund.lasso$Abundance)
summary(fit)
dev.off()

####ok and if you want to make sure the connection is EVER present:
#lassoplots connection was ever present
#for each cluster-patient-speciesgroup combo from above...can you get the abundances at each timepoint for that comb
onlyspecies.arg <- subset(df.arg.agg.species, Speciesname!="s__.")
#bring in somehow all the timepoints before you aggregate
#ugh this is missing the zeroes for timepoints
#first merge abundance data with a file that has all timepoints vs all species
samples = unique(abund.species.agg$Sample)
species = unique(abund.species.agg$Species)
allcombinations = expand.grid(samples,species)
colnames(allcombinations)<-c("Sample","Species")
abund.species.agg.allcombos <- merge(allcombinations, abund.species.agg, by=c("Sample","Species"),all.x=T)
abund.species.agg.allcombos$Patient<-as.vector(sapply(strsplit(as.character(abund.species.agg.allcombos$Sample),"-"), function(x) x[[1]]))
abund.species.agg.allcombos$Abundance[is.na(abund.species.agg.allcombos$Abundance)]<-0
arg.speciesabundance<-merge(onlyspecies.arg, abund.species.agg.allcombos, by=c("Patient","Species"),all.x=T)

#bring in the ARG abundances based on the sample
argrpkm.speciesabundance <- merge(arg.speciesabundance, df5.melt, by.x=c("Sample","Patient","Cluster"),by.y=c("Sample","Patient","Cluster"))
argrpkm.speciesabund.lasso<-aggregate(Abundance~Sample +Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Sample+Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df <- merge(argrpkm.speciesabund.lasso, toparg_clusters, by="Cluster", all.x=T)
df <-merge(known_cat,df,by="Mergedname")
df.to.plot <- subset(argrpkm.speciesabund.lasso)
df.to.plot<- df
pdf("TopARGS_vs_Species_abundances_lasso_samples-presentever.pdf",height=7,width=9.5,useDingbats=F)
plot22 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_text(aes(label=Trim_name, cex=2),hjust=-.1, vjust=0)+
        geom_point(aes(col=Patient, size = Count))+
        xlab("Contributing Species Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=pats, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15,20))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot22)
fit = lm(df.to.plot$RPKM~df.to.plot$Abundance)
summary(fit)
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="pearson")
cor.test(df.to.plot$RPKM,df.to.plot$Abundance,method="spearman")
dev.off()


#what if you want to keep the samples separate
sub_argcount <- top10
sub_argcount$connectsum <-NULL
sub_argcount$Trim_name<-NULL
sub_argcount$Confidence <-NULL
sub_argcount$ARG_name<-NULL
sub_argcount$Specific_name<-NULL
sub_argcount$Top_ARG<-NULL
#keep the sample in
df_melt2 <- melt(sub_argcount, id = c("Sample","patients", "Cluster", "Mergedname","Category"))
df_melt.pos <-subset(df_melt2, value>0)
df_melt.pos$levels <-as.vector(sapply(strsplit(as.character(df_melt.pos$variable),"[.][.]"), function(x) list(trimws(x))))
df_melt.pos$species <- sapply(df_melt.pos$levels, function(x) paste(x[1:7], collapse=";"))
df_melt.pos$speciesname <- sapply(df_melt.pos$levels, function(x) paste(x[7:7], collapse=";"))
#species
#aggregate keeping clusters, patients, and taxa
df.arg.agg.species <-aggregate(value ~Sample+patients+Cluster+Mergedname+Category+species+speciesname, data =  df_melt.pos, sum)
df.arg.agg.species$value <- (df.arg.agg.species$value>0) + 0
colnames(df.arg.agg.species)<-c("Sample","Patient","Cluster","Mergedname","Category","Species", "Speciesname","Count")
df.arg.agg.species.reduced<-aggregate(Count~Sample+Patient+Species+Speciesname,data=df.arg.agg.species,sum)
colnames(df.arg.agg.species.reduced)<-c("Sample","Patient","Species","Speciesname","Count")

pdf("TopARGS_vs_Species_abundances_samples.pdf",useDingbats=F)
species.argrpkm <- merge(df.arg.agg.species, df5.melt, by=c("Sample","Patient","Cluster"))
argrpkm.speciesabundance<-merge(species.argrpkm, abund.species.agg, by=c("Sample","Patient","Species"))
plot23 <-ggplot(argrpkm.speciesabundance, aes(x = Abundance, y = RPKM))+
        geom_point(aes(col=Patient))+
        xlab("Taxa Abundance")+
        ylab("ARG Abundance (RPKM)")+
        scale_color_manual(label=patients, values=patient.col)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot23)
dev.off()

#next goal is to take these top args and add all of the taxa that represent them together...so that you only have one line for each arg...witha summed species abundance
#is it a straight line? then the species are representative of total populations
#is it not a straight line?

pdf("TopARGS_vs_Species_abundances_lasso_samples-present.pdf",height=7, width=9.5,useDingbats=F)
species.argrpkm <- merge(df.arg.agg.species, df5.melt, by=c("Sample","Patient","Cluster"))
onlyspecies.argrpkm <- subset(species.argrpkm, Speciesname!="s__.")
argrpkm.speciesabundance<-merge(onlyspecies.argrpkm, abund.species.agg, by=c("Sample","Patient","Species"))
argrpkm.speciesabund.lasso<-aggregate(Abundance~Sample+Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.count<-aggregate(Count~Sample+Patient+Cluster+RPKM,data=argrpkm.speciesabundance,sum)
argrpkm.speciesabund.lasso$Count <- argrpkm.speciesabund.count$Count
df <- merge(argrpkm.speciesabund.lasso, toparg_clusters, by="Cluster", all.x=T)
df <-merge(known_cat,df,by="Mergedname")
df.to.plot<- df
plot24 <-ggplot(df.to.plot, aes(x = Abundance, y = RPKM))+
        geom_text(aes(label=Trim_name, cex=2),hjust=-.1, vjust=0)+
        geom_point(aes(col=Category, size = Count))+
        xlab("Contributing Species Abundance (Sum(Median(Species Contig RPKMs)))")+
        ylab("ARG Abundance (RPKM)")+
        #scale_color_manual(label=, values=patient.col)+
        scale_size(name= "Number of Species\nContributing",breaks = c(1,5,10,15))+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(plot24)
print("this is the correlation information for the TopARGS_vs_Species_abundances_lasso_samples-present:")
fit <-lm(argrpkm.speciesabund.lasso$RPKM~argrpkm.speciesabund.lasso$Abundance)
summary(fit)
cor.test(argrpkm.speciesabund.lasso$RPKM,argrpkm.speciesabund.lasso$Abundance)
cor.test(argrpkm.speciesabund.lasso$RPKM,argrpkm.speciesabund.lasso$Abundance,method="spearman")
dev.off()


subset(argrpkm.speciesabund.lasso, RPKM>100)

