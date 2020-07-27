library(ggplot2)
library(scales)

setwd("/workdir/users/agk85/CDC2/bins")

metadata = read.csv("/workdir/users/agk85/CDC2/bins/all_das_checkm_kraken.txt", header = F,sep='\t')
colnames(metadata)<-c("Bin","Checkm_taxa","Completeness","Contamination","Heterogeneity","Genome_size","Kraken_completeness","Taxid","Rank_name","Rank","Kraken_contamination","Taxonomy")

#colors = checkm completeness
#shape = rank_name

pdf("Contamination_comparison.pdf", height = 10, width = 10)
ggplot(metadata, aes(Kraken_contamination, Contamination))+
geom_point(aes(shape=Rank_name, color=Completeness))+
scale_x_sqrt(labels=comma)+
scale_y_sqrt(labels=comma)+
scale_shape_manual(values = c(0,1,2,3,4,5,6,8))+
scale_color_gradient(low = "#132B43", high = "#56B1F7",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "color")
dev.off()

pdf("Contamination_comparison2.pdf", height = 10, width = 10) 
ggplot(metadata, aes(Kraken_contamination, Contamination))+
geom_point(aes(fill=Rank), shape=21)+
scale_x_sqrt(labels=comma)+
scale_y_sqrt(labels=comma)+
#scale_shape_manual(values = c(0,1,2,3,4,5,6,8))+
scale_fill_gradient(low = "white", high = "maroon",
  space = "Lab", na.value = "grey50", guide = "colourbar")
dev.off()


dev.off()











