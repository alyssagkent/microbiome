library(ggplot2)
library(RColorBrewer)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/histograms")


#color generation
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 7
kewl<-col_vector[1:(n)]
pie(rep(1,n), col=kewl)

#import metaphlan
metaphlan <- read.table("/workdir/users/agk85/CDC/metaphlan/combo3/temp/Mgm_metaphlan.txt",sep="\t",header=F,row.names=1)
metaphlan<- as.matrix(metaphlan)
n <- as.vector(metaphlan[1,1:ncol(metaphlan)])
colnames(metaphlan) <- n
enterobacteriaceae <- as.numeric(metaphlan["k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacteriales|f__Enterobacteriaceae",])

#import ARG RPKM
infile2=read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95.txt",header =T,row.names=1)
total_rpkm <- rowSums(infile2)
tot_rpkm <- total_rpkm[order(names(total_rpkm))]

argentero <- data.frame(names(tot_rpkm), tot_rpkm, enterobacteriaceae)
colnames(argentero)<-c("Sample","ARGrpkm","Enterobacteriaceae")
argentero$Patient <- as.factor(as.vector(sapply(strsplit(as.character(argentero$Sample),"-"), function(x) x[[1]])))

pdf("ARG_enterobacteriaceae.pdf",height=5, width=7,useDingbats=F)
ggplot(data = argentero)+
xlab("Relative Abundance of Enterobacteriaceae")+
ylab("Total ARG RPKM")+
geom_point(aes(x=Enterobacteriaceae, y=ARGrpkm,col=Patient))+
scale_colour_manual(values=kewl)
dev.off()
