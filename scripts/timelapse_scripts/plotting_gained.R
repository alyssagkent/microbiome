#this program will plot the histograms for one sample
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(colorspace)
library(reshape2)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC2/bins/timelapse")
args <- commandArgs(trailingOnly = TRUE)
thresh=args[1]
genetype = args[2] #arg or mge
taxlevel = args[3]

inhandle = paste('timelapse_', genetype, '_org_',taxlevel,'_',thresh,'_taxa_species_gained.txt',sep='')
df = read.csv(inhandle,header=T, sep="\t")

pdf(paste("Gained_taxa_",genetype, '_org_', taxlevel, '_',thresh,'.pdf',sep=''),height=30,width=20)
ggplot(df, aes(x=Taxa,fill=T2))+
 geom_bar(stat="count",color="black")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 facet_grid(T1~., scales = "free_y")
dev.off()


argfile = read.csv("/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech.txt",header=T, sep="\t")

mergedf = merge(df, argfile, by.x="Gene", by.y="Cluster")

pdf(paste("Gained_genes_type_",genetype, '_org_', taxlevel, '_',thresh,'.pdf',sep=''),height=30,width=20)
ggplot(mergedf, aes(x=Mechanism,fill=T2))+
 geom_bar(stat="count",color="black")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 facet_grid(T1~., scales = "free_y")
dev.off()

pdf(paste("Gained_genes_name_",genetype, '_org_', taxlevel, '_',thresh,'.pdf',sep=''),height=30,width=20)
ggplot(mergedf, aes(x=Name,fill=T2))+
 geom_bar(stat="count",color="black")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 facet_grid(T1~., scales = "free_y")
dev.off()




pdf(paste("Gained_genes_card_",genetype, '_org_', taxlevel, '_',thresh,'.pdf',sep=''),height=30,width=20)
ggplot(mergedf, aes(x=CARD,fill=T2))+
 geom_bar(stat="count",color="black")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 facet_grid(T1~., scales = "free_y")
dev.off()
