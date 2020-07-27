library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(scales)
library(circlize)
library(ggplot2)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())

setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/abundances")

drug <- read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/ARG_ofinterest_annotations.txt",sep="\t",header=T)



#go to taxonomic categories
df <- merge(df_melt2, drug, by=c("Cluster"))
df2 <- subset(df, value>0)
df_levofloxacin <- subset(df2, Drug=="levofloxacin")
df_piptazo <- subset(df2, Drug=="piperacillin-tazobactam" | Drug == "piperacillin")
df_tmpsmx <- subset(df2, Drug=="trimethoprim" | Drug == "Sulfamethoxazole")
df_levofloxacin$



pdf("Levofloxacin_gene_distribution.pdf",height=15, width=5)
ggplot(df_levofloxacin, aes(variable, fill=patients))+
geom_histogram(position = "stack", stat="count")+
scale_fill_manual(values=kewl)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf("Pip-Tazo_gene_distribution.pdf",height=15, width=5)
ggplot(df_piptazo, aes(variable, fill=patients))+
geom_histogram(position = "stack", stat="count")+
scale_fill_manual(values=kewl)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

pdf("TMP-SMX_gene_distribution.pdf",height=15, width=10)
ggplot(df_tmpsmx, aes(variable, fill=patients))+
geom_histogram(position = "stack", stat="count")+
scale_fill_manual(values=kewl)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()





library("dplyr")
data %>% group_by(Variable) %>% summarize(count=n())

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

df6 <- merge(df5.melt, drug, by=c("Cluster"))

scale_color_manual(values=kewl)
dev.off()

pdf("Lachnospiraceae_metaphlan.pdf")
ggplot(data=lachno_df, aes(TP, Taxa, group=Patient))+
geom_line(aes(col=Patient))+
scale_color_manual(values=kewl)
dev.off()

tp_dp <- read.csv("TP_DP.txt",sep="\t",header=T)
drugs <- read.csv("drug_data.txt",sep="\t",header=T)
drugs$hicd <- drugs$last_metagenome_dp - drugs$Start

drug_subset <- subset(drugs, Study.Patient %in% unique(Patient) & hicd>0)
drug_sub <- merge(drug_subset, num, by.y=c("Patient"), by.x=c("Study.Patient"))

num = data.frame(c("B314","B316","B320","B331","B335","B357","B370"),c(7,6,5,4,3,2,1))
colnames(num)<-c("Patient","Patient_Number")


drugbar <- ggplot()+
geom_rect(data=num, mapping=aes(xmin=-5, xmax =60, ymin = Patient_Number-0.5,ymax=Patient_Number+0.5, fill=Patient),alpha=0.2)+
scale_fill_manual(values=c(kewl,c("#B6D37C","#80A88A","#2F6A39","#928D8D")))+
ylab("Patient")+
geom_rect(data=drug_sub, mapping=aes(xmin=Start, xmax=Stop, ymin=Patient_Number-0.5+(Antibacterial_number-1)*1/3, ymax=Patient_Number-0.5 + (Antibacterial_number*1/3), fill=Antibacterial.1), color="black", alpha=0.5)+
theme(legend.position="none")



top20families <-c("k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae",
"k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae",
"k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae",
"k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae",
"k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae",
"k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Porphyromonadaceae",
"k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selenomonadales|f__Veillonellaceae",
"k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Caudovirales|f__Siphoviridae",
"k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Streptococcaceae",
"k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Lactobacillaceae",
"k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae",
"k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Enterococcaceae",
"k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Coriobacteriales|f__Coriobacteriaceae",
"k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Pasteurellales|f__Pasteurellaceae",
"k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Rikenellaceae",
"k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Clostridiaceae",
"k__Bacteria|p__Firmicutes|c__Erysipelotrichia|o__Erysipelotrichales|f__Erysipelotrichaceae",
"k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae",
"k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Verrucomicrobiaceae",
"k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selenomonadales|f__Acidaminococcaceae")
figs = vector("list",20)
for (i in seq(length(top20families))){
family = top20families[i]
print(family)
print(i)
df1 <- data.frame(Sample, Patient, TP, as.numeric(metaphlan[family,]))
df <- merge(df1, tp_dp, by=c("Patient","TP"))
fname = strsplit(family, "[|]")[[1]][[5]]
colnames(df)<- c("Patient","TP","Sample","Taxa","DP")
pdf(paste("DPA_abundance_", fname,"_metaphlan.pdf",sep=""))
figure = ggplot(data=df, aes(DP, Taxa, group=Patient))+
ylab(paste("Relative Abundance of ", fname, sep=""))+
xlab("Days Post Admission")+
xlim(c(min(df$DP), 60))+
geom_line(aes(col=Patient))+
scale_color_manual(values=kewl)+
theme(legend.position="none")
print(figure)
dev.off()
figs[[i]] = figure
}

pdf(paste("All_top20_families_abundance_metaphlan_dpa.pdf",sep=""), height = 20, width = 30)
grid.arrange(drugbar, drugbar, drugbar, drugbar, drugbar, 
	figs[[1]], figs[[2]], 
        figs[[3]], figs[[4]],
        figs[[5]], figs[[6]],
        figs[[7]], figs[[8]],
        figs[[9]], figs[[10]],
        figs[[11]], figs[[12]],
        figs[[13]], figs[[14]],
        figs[[15]], figs[[16]],
	figs[[17]], figs[[18]],
	figs[[19]], figs[[20]],
        ncol=5, nrow=5, widths=c(rep(4,5)), heights=c(2,rep(5,4)))
dev.off()



figs = vector("list",20)
for (i in seq(length(top20families))){
family = top20families[i]
print(family)
print(i)
dfprime <- data.frame(Sample, Patient, TP, as.numeric(metaphlan[family,]))
df <- merge(dfprime, tp_dp, by=c("Patient","TP"))
fname = strsplit(family, "[|]")[[1]][[5]]
colnames(df)<- c("Patient","TP","Sample","Taxa","DP")
pdf(paste("TP_abundance_", fname,"_metaphlan.pdf",sep=""))
figure = ggplot(data=df, aes(TP, Taxa, group=Patient))+
ylab(paste("Relative Abundance of ", fname, sep=""))+
xlab("Timepoints")+
geom_line(aes(col=Patient))+
scale_color_manual(values=kewl)
print(figure)
dev.off()
figs[[i]] = figure
}

geom_segment(data=drugset,aes(x=DP.start, y=0, xend=DP.stop, yend=0), colour="purple", alpha=.2, size=100)



pdf(paste("All_top20_families_abundance_metaphlan_tp.pdf",sep=""), height = 20, width = 30)
grid.arrange(figs[[1]], figs[[2]], 
        figs[[3]], figs[[4]],
        figs[[5]], figs[[6]],
        figs[[7]], figs[[8]],
        figs[[9]], figs[[10]],
        figs[[11]], figs[[12]],
        figs[[13]], figs[[14]],
        figs[[15]], figs[[16]],
        figs[[17]], figs[[18]],
        figs[[19]], figs[[20]],
        ncol=5, nrow=4, widths=c(rep(4,5)), heights=c(rep(5,4)))
dev.off()




#so what if you use the bacterial average family abundances you used for the other figures
#import max abundance
abund <-read.table("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/all_taxonomic_abundances_basehic_nophage.txt",sep='\t', header = F)
colnames(abund)<-c("Scaffold","Sample","Patient","Abundance","Organism")
#aggregate to total abundance of each organism at family level
abund$Levels <-as.vector(sapply(strsplit(as.character(abund$Organism),"; "), function(x) list(trimws(x))))
abund$Family <- sapply(abund$Levels, function(x) paste(x[5:5], collapse=";"))
abund.fam.agg <-aggregate(Abundance ~ Sample +  Patient + Family, data = abund, sum)
abund.tot.enteros <- subset(abund.fam.agg, Family == "f__Enterobacteriaceae")


abund$Species <- sapply(abund$Levels, function(x) paste(x[7:7], collapse=";"))
abund_onlyspecies <- subset(abund, Species!="s__;")
abund.species.agg <-aggregate(Abundance ~ Sample + Patient + Species + Family, data = abund_onlyspecies, median)
abund.fam.avg.agg <-aggregate(Abundance ~ Sample + Patient + Family, data=abund.species.agg, sum) 

#aggregate to all 
abund.all.agg <-aggregate(Abundance ~Sample + Patient, data = abund.fam.agg, sum)
#merge with the metaphlan to assess the correlation
merge_df1<- merge(metaphlan_entero, abund.tot.enteros, by="Sample",all.x=T)
merge_df1$Abundance[is.na(merge_df1$Abundance)]<-0
merge_df1$All_Abundance <- abund.all.agg$Abundance
merge_df1$Ratio <- merge_df1$Abundance/merge_df1$All_Abundance

merge_df <- merge(merge_df1, abund.avg.enteros, by="Sample", all.x=T)
merge_df$Average_Abundance[is.na(merge_df$Average_Abundance)]<-0

plot(merge_df$Metaphlan_Enterobacteriaceae, merge_df$Abundance)
plot(merge_df$Metaphlan_Enterobacteriaceae, merge_df$Ratio)
plot(merge_df$Metaphlan_Enterobacteriaceae, merge_df$Average_Abundance)
#import ARG RPKM
infile2=read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95.txt",header =T,row.names=1)
total_rpkm <- rowSums(infile2)
tot_rpkm <- total_rpkm[order(names(total_rpkm))]

argentero <- data.frame(names(tot_rpkm), tot_rpkm, 100*merge_df$Ratio, merge_df$Abundance, merge_df$Average_Abundance,merge_df$Metaphlan_Enterobacteriaceae)
colnames(argentero)<-c("Sample","ARGrpkm","Relative_Enterobacteriaceae","Enterobacteriaceae","Average_Enterobacteriaceae","Metaphlan_Enterobacteriaceae")
argentero$Patient <- as.factor(as.vector(sapply(strsplit(as.character(argentero$Sample),"-"), function(x) x[[1]])))



#############################################
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
abund$Family <- sapply(abund$Levels, function(x) paste(x[1:5], collapse=";"))
#average all of the contigs based on their sample/patient/species
abund.species.agg <-aggregate(Abundance ~ Sample +  Patient + Species, data = abund, median)
#take the maximal abundance of a particular sample
abund.species.agg.patient <-aggregate(Abundance ~ Patient + Species, data =abund.species.agg, max)


abund.agg.species <- aggregate(Abundance ~ Sample +  Patient + Species+Family, data = abund, median)
abund.agg.family <- aggregate(Abundance ~ Sample + Patient + Family, data = abund.agg.species, sum)
abund.agg.family$TP <- as.vector(sapply(strsplit(as.character(abund.agg.family$Sample),"-"), function(x) x[[2]])

ggplot(data=abund.agg.family, aes(TP, Abundance, group=Family))+
geom_line()+
facet_grid(Patient~., scale="free_y")

sub = subset(abund.agg.family, Family =="k__Bacteria;p__Firmicutes;c__Erysipelotrichia;o__Erysipelotrichales;f__Erysipelotrichaceae")
ggplot(data=sub, aes(TP, Abundance, group=Family))+
geom_line(aes(col=Family))+
facet_grid(Patient~., scale="free_y")
"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Shigella;s__Shigella_flexneri."
"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Shigella;s__Shigella_sp._PAMC_28760."

sub = subset(abund.agg.family, Species =="k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Shigella;s__Shigella_dysenteriae.")
ggplot(data=sub, aes(TP, Abundance, group=Family))+
geom_line(aes(col=Family))+
facet_grid(Patient~., scale="free_y")


