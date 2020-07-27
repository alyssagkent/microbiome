library(ggplot2)
library(RColorBrewer)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/taxonomic_distributions")


#color generation
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 7
kewl<-col_vector[1:(n)]
#pie(rep(1,n), col=kewl)

#import metaphlan
metaphlan <- read.table("/workdir/users/agk85/CDC/metaphlan/combo3/temp/Mgm_metaphlan.txt",sep="\t",header=F,row.names=1)
metaphlan<- as.matrix(metaphlan)
n <- as.vector(metaphlan[1,1:ncol(metaphlan)])
colnames(metaphlan) <- n
enterobacteriaceae <- as.numeric(metaphlan["k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae",])
metaphlan_entero <- data.frame(colnames(metaphlan), enterobacteriaceae)
colnames(metaphlan_entero)<- c("Sample","Metaphlan_Enterobacteriaceae")

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
abund.avg.enteros <- subset(abund.fam.avg.agg, Family == "f__Enterobacteriaceae")
colnames(abund.avg.enteros)<-c("Sample","Patient","Family","Average_Abundance")
abund.avg.enteros$Patient<-NULL
abund.avg.enteros$Family<-NULL

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

pdf("ARG_enterobacteriaceae_summedaverage_rpkm.pdf",height=5, width=7,useDingbats=F)
ggplot(data = argentero)+
xlab("Abundance of Enterobacteriaceae (Sum(Median(Species RPKM)))")+
ylab("Total Abundance of Antibiotic Resistance Genes (RPKM)")+
geom_point(aes(x=Average_Enterobacteriaceae, y=ARGrpkm,col=Patient))+
scale_colour_manual(values=kewl)
dev.off()

corstat1 <- cor.test(argentero$Average_Enterobacteriaceae, argentero$ARGrpkm)


#
#linear mixed effect modeling #alyssa's 1st attempt
library(lme4)
library(lmerTest)

library(nlme)
m <- lme(ARGrpkm ~ 1+ Average_Enterobacteriaceae, random= ~ 1 + Average_Enterobacteriaceae|Patient, data=argentero)
summary(m)

ggplot(argentero, aes(Average_Enterobacteriaceae, ARGrpkm)) +
  geom_point(color="grey") + 
  geom_smooth(method="lm", se = FALSE) +
  geom_line(aes(y=predict(m)), color="red") +
  facet_wrap( ~ Patient) 




m1 <- lmer(ARGrpkm ~ Relative_Enterobacteriaceae + (1+Relative_Enterobacteriaceae|Patient), argentero)
summary(m1)
m2 <- lmer(ARGrpkm ~ Average_Enterobacteriaceae + (1+Average_Enterobacteriaceae|Patient), argentero)
summary(m2)
m3 <- lmer(ARGrpkm ~ Metaphlan_Enterobacteriaceae + (1+Metaphlan_Enterobacteriaceae|Patient), argentero)
summary(m3)

m4 <- lmer(ARGrpkm ~ Relative_Enterobacteriaceae + (1|Patient), argentero)
summary(m4)
m5 <- lmer(ARGrpkm ~ Average_Enterobacteriaceae + (1|Patient), argentero)
summary(m5)
m6 <- lmer(ARGrpkm ~ Metaphlan_Enterobacteriaceae + (1|Patient), argentero)
summary(m6)





pdf("ARG_metaphlan_enterobacteriaceae.pdf",height=5, width=7,useDingbats=F)
ggplot(data = argentero)+
xlab("Relative Abundance of Enterobacteriaceae")+
ylab("Total Abundance of Antibiotic Resistance Genes (RPKM)")+
geom_point(aes(x=Metaphlan_Enterobacteriaceae, y=ARGrpkm,col=Patient))+
scale_colour_manual(values=kewl)
dev.off()











pdf("ARG_enterobacteriaceae_relative_rpkm.pdf",height=5, width=7,useDingbats=F)
ggplot(data = argentero)+
xlab("Relative Abundance of Enterobacteriaceae")+
ylab("Total Abundance of Antibiotic Resistance Genes (RPKM)")+
geom_point(aes(x=Relative_Enterobacteriaceae, y=ARGrpkm,col=Patient))+
scale_colour_manual(values=kewl)
dev.off()
corstat2 <- cor.test(argentero$Relative_Enterobacteriaceae, argentero$ARGrpkm)

stat1 <- c(corstat1$estimate, corstat1$p.value)
stat2 <- c(corstat2$estimate, corstat2$p.value)
r = rbind(stat1, stat2)
write.csv(r, "arg_enterobacteriaceae_statistics.txt", quote=F)


#pdf("ARG_enterobacteriaceae_total_rpkm.pdf",height=5, width=7,useDingbats=F)
#ggplot(data = argentero)+
#xlab("Total Abundance of Enterobacteriaceae (RPKM)")+
#ylab("Total Antibiotic Resistance Gene (RPKM)")+
#geom_point(aes(x=Enterobacteriaceae, y=ARGrpkm,col=Patient))+
#scale_colour_manual(values=kewl)
#dev.off()
