
#make sure you get the header off of the file after catting
#cat *_taxonomy_overlap.txt | sort | uniq > Taxonomy_overlap.txt
#grep -v '^beginning name' | sort | uniq > Taxonomy_overlap.txt (might get you there!)
#vim Taxonomy_overlap.txt, scroll and remove header

#makes a figure of the Taxonomy overlap between gaemr, kraken, original bestorg ----only works if 
library(ggplot2)
library(reshape2)
set.seed(42)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())

setwd("/workdir/users/agk85/CDC2/combo_tables/metagenomes/")
df <-read.csv("Taxonomy_overlap.txt",header=F, sep="\t")
colnames(df)<-c("Sample","Group","Kingdom","Phylum","Class","Order","Family","Genus","Species","Total_contigs","Kingdom_proportion","Phylum_proportion","Class_proportion","Order_proportion","Family_proportion","Genus_proportion","Species_proportion")

dfsub <- subset(df, select=c("Sample","Group","Kingdom_proportion","Phylum_proportion","Class_proportion","Order_proportion","Family_proportion","Genus_proportion","Species_proportion"))

goodgroups <- c("gb","kb","kg","kgb")
dfsubset <- subset(dfsub, Group %in% goodgroups)

dfmelt <- melt(dfsubset)
colnames(dfmelt) <- c("Sample","Group","Taxagroup","Proportion")
pdf("Taxonomy_overlap.pdf",height=15,width=30)
ggplot(dfmelt,aes(x=interaction(Group,Taxagroup),y= Proportion))+
geom_boxplot(aes(fill=Group),outlier.shape=NA)+
facet_grid(. ~ Taxagroup,scales="free_x")+
geom_line(aes(group=interaction(Sample,Taxagroup)),
	alpha = 0.5, colour ="black")
#geom_point()+
#geom_jitter(width = 0.25)
dev.off()


pdf("Taxonomy_overlap.pdf",height=5,width=15)
ggplot(dfmelt,aes(x=Group,y= Proportion))+
geom_boxplot(aes(fill=Group),outlier.shape=NA)+
facet_grid(. ~ Taxagroup,scales="free_x")+
scale_fill_manual(values=c("gb"="darkorange1","kb"="darkcyan","kg"="darkslateblue","kgb"="darkorange4"),
breaks=c("gb","kb", "kg", "kgb"),
                       labels=c("Gaemr-Original","Kraken-Original", "Kraken-Gaemr", "Kraken-Gaemr-Original"))+
geom_point(size=.5,position=position_jitter(.25))
dev.off()


