#withinbetween smillie taxonomic levels
library(ggplot2)
library(reshape2)
library(ggpubr)

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC2/bins/hgt_comparisons")

args <- commandArgs(trailingOnly = TRUE)
thresh=args[1] #2 or 5
genetype = args[2] #arg or mge

indata = read.table(paste("hgt_comparisons_alllevels_allorgs_", genetype, "_", thresh, ".txt",sep=""),header=F, sep="\t")
colnames(indata)<-c("Patient1","Patient2","Linked","Total","Withinbetween","Taxalevel")
indata$Withinbetween <- factor(indata$Withinbetween, levels = c("within","between"))
indata$Taxalevel <-factor(indata$Taxalevel, levels=c("genus","family","order","class","phylum","kingdom"))
indata$Ratio <- indata$Linked/indata$Total
#test before you change the TaxaLevel names

g.test = wilcox.test(Ratio~Withinbetween, data=subset(indata, Taxalevel =="genus"), exact=F)
f.test = wilcox.test(Ratio~Withinbetween, data=subset(indata, Taxalevel =="family"),exact=F)
o.test = wilcox.test(Ratio~Withinbetween, data=subset(indata, Taxalevel =="order"),exact=F)
c.test = wilcox.test(Ratio~Withinbetween, data=subset(indata, Taxalevel =="class"),exact=F)
p.test = wilcox.test(Ratio~Withinbetween, data=subset(indata, Taxalevel =="phylum"),exact=F)
k.test = wilcox.test(Ratio~Withinbetween, data=subset(indata, Taxalevel =="kingdom"),exact=F) 
g.test
f.test
o.test
c.test
p.test
k.test



levels(indata$Taxalevel)<-c("Between\nSpecies","Between\nGenera","Between\nFamilies","Between\nOrders","Between\nClasses","Between\nPhyla")
cm = compare_means(Ratio ~ Withinbetween, data = indata, method = "wilcox.test", p.adjust.method="bonferroni",paired = FALSE,group.by = "Taxalevel")
cm
write.table(cm, paste(genetype, "_",thresh,"_hgt_statistics.txt",sep=""), sep="\t")

pdf(paste(genetype, "_",thresh,"_hgt_comparisons_all.pdf",sep=""), useDingbats=F,height = 5, width =8)
ggplot(data=indata, aes(x=Withinbetween,y=Ratio*100))+
	geom_boxplot(aes(fill=Withinbetween),outlier.shape = NA)+
	ylab("HGT rate (per 100 comparisons)\nof genes")+
	xlab("Taxonomic Levels")+
	scale_fill_manual(values=c(title="","within"="red","between"="blue"), labels=c("Within Patient","Between Patients"))+
	facet_grid(.~Taxalevel, scales="free")+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.1))+
	geom_jitter(position=position_jitter(0.2))
dev.off()

write.csv(indata, paste(genetype, "_",thresh,"_hgt_comparisons_alllevels_ratio.txt",sep=""))

