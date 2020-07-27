library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(igraph)
set.seed(42)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())
args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
contacts = args[2]
genetype = args[3]
line_scaling = as.numeric(args[4])

#infile = '/workdir/users/agk85/CDC2/bins/arg_das_5_taxa_taxa_counts.txt'
#contacts = 5
###########
#import data
###########
setwd("/workdir/users/agk85/CDC2/bins")
df <-read.table(infile,header=T, row.names=1,sep="\t")
# reduce number of counts by 2 because of the doubling
df$Genes <- df$Genecount/2

################
#colors
################
arg.palette <-colorRampPalette(brewer.pal(12,"Set3"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 9
kewl<-col_vector[60:(n+60)]
################################
mypalette  <-colorRampPalette(c("red", "blue"))(n=50)
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}



df2 <- subset(df, Level =="s__" & Patient != 'all', select =c(Patient, Taxa1, Taxa2, Genes))
df2$Connections = paste(df2$Taxa1, df2$Taxa2)
df3 <- subset(df2, select = c(Patient, Connections, Genes))
df4 <- dcast(df3, Patient ~ Connections)
#replace NAs with 0's
df4[is.na(df4)]<-0
rownames(df4) <- df4$Patient
df4$Patient <- NULL
df5 <- data.matrix(df4)
df5.dist <- vegdist(as.matrix(df5), method="jaccard", binary=T)
df5.dist.mat <- as.matrix(df5.dist)
write.table(df5.dist.mat, paste("patient_jaccard_distance_connections_",contacts,"_",genetype,".txt",sep=""), sep="\t", quote=FALSE)


df5.dist.mat[df5.dist.mat==0]<-NA
pdf(paste("patient_jaccard_distance_connections_",contacts,"_",genetype,".pdf",sep=""), height=20, width=20,useDingbats=FALSE)
heatmap.2((df5.dist.mat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
        dendrogram="both", Rowv=T,Colv=T,key=TRUE,symbreaks=FALSE,
	density.info="none", trace="none", na.color = "white",
	labCol=colnames(df5.dist.mat), labRow=rownames(df5.dist.mat),
	cexCol=2, col=mypalette,keysize=2,key.par = list(cex=2))
dev.off()
#bray_curtis#########################

df5.dist <- vegdist(as.matrix(df5), method="bray")
df5.dist.mat <- as.matrix(df5.dist)
write.table(df5.dist.mat, paste("patient_bray_distance_connections_",contacts,"_",genetype,".txt",sep=""), sep="\t", quote=FALSE)

df5.dist.mat[df5.dist.mat==0]<-NA
pdf(paste("patient_bray_distance_connections_",contacts,"_",genetype,".pdf",sep=""), height=20, width=20,useDingbats=FALSE)
heatmap.2((df5.dist.mat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
        dendrogram="both", Rowv=T,Colv=T,key=TRUE,symbreaks=FALSE,
        density.info="none", trace="none",na.color = "white",
        labCol=colnames(df5.dist.mat), labRow=rownames(df5.dist.mat),
        cexCol=2, col=mypalette,keysize=2,key.par = list(cex=2))
dev.off()

