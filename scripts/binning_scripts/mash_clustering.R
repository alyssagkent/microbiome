library(ggplot2)
library(gplots)

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())

setwd("/workdir/users/agk85/CDC2/bins/")
d = read.csv("/workdir/users/agk85/CDC2/bins/hq_bins_dist.txt", stringsAsFactors=FALSE, header= F, row.names=1, sep='\t')
n = as.character(unlist(d[1,]))

colnames(d)<-n
d<-d[-1,]

colnames(d)=gsub(".fa","",colnames(d))
rownames(d)=gsub(".fa","",rownames(d))

dmat <- data.matrix(d)

metadata = read.csv("/workdir/users/agk85/CDC2/bins/all_das_checkm_kraken.txt", header = F,sep='\t')
colnames(metadata)<-c("Bin","Checkm_taxa","Completeness","Contamination","Heterogeneity","Genome_size","Kraken_completeness","Taxid","Rank_name","Rank","Kraken_contamination","Taxonomy")



df <- merge(metadata, dmat, by.x="Bin",by.y="row.names")
df.metadata <- df[,1:12]
df.data <- df[,13:ncol(df)]
rownames(df.metadata)<-df$Bin
rownames(df.data)<-df$Bin
df.datamat <- data.matrix(df.data)



pdf("/workdir/users/agk85/CDC2/bins/hq_bins_heatmap_mash.pdf", height=100, width=100)
heatmap.2(df.datamat, labRow=paste(df.metadata$Taxonomy, " ", df.metadata$Bin,sep=""), labCol=df.metadata$Taxonomy, 
trace="none", cexRow=.3,cexCol=0.3,lhei=c(1, 10),lwid=c(1,10),margins=c(40,40))
dev.off()


df.metadata.species <- subset(df.metadata, Rank_name=="species")
speciesbins <- df.metadata.species$Bin
df.data.species <- df.datamat[rownames(df.datamat) %in% speciesbins,colnames(df.datamat) %in% speciesbins ]


pdf("/workdir/users/agk85/CDC2/bins/hq_bins_species_heatmap_mash.pdf", height=100, width=100)
heatmap.2(df.data.species, labRow=paste(df.metadata.species$Taxonomy, " ", df.metadata.species$Bin,sep=""),
labCol=paste(df.metadata.species$Taxonomy, " ", df.metadata.species$Bin,sep=""), 
trace="none", cexRow=.8,cexCol=0.8,lhei=c(.5, 10),lwid=c(.5,10),margins=c(60,60))
dev.off()


df.metadata.krakencont <- subset(df.metadata, Kraken_contamination<10)
krakencontbins <- df.metadata.krakencont$Bin
df.data.krakencont <- df.datamat[rownames(df.datamat) %in% krakencontbins,colnames(df.datamat) %in% krakencontbins]
labels = paste(df.metadata.krakencont$Taxonomy, " ", df.metadata.krakencont$Bin,sep="")

pdf("/workdir/users/agk85/CDC2/bins/hq_bins_krakencont_heatmap_mash.pdf", height=100, width=100)
heatmap.2(df.data.krakencont, labRow=labels,labCol=labels,
trace="none", cexRow=.8,cexCol=0.8,lhei=c(.5, 10),lwid=c(.5,10),margins=c(60,60))
dev.off()
###############################nmds
library(MASS)
library(ggrepel)
# Read in the distance file, aaf.dist for example
distance <-as.dist(df.datamat)
# Turn it into a matrix
colnames(distance) <- rownames(distance)
dis_matrix <- as.matrix(distance)
# Collaps into two dimension using NMDS argorithm
dis_nmds <- isoMDS(dis_matrix)
# Plot 
points <- as.data.frame(dis_nmds$points)
colnames(points) <- c('NMDS1','NMDS2')


#combine with metadata
nmds.df <- merge(points, metadata, by.x="row.names",by.y="Bin")
#get taxa names
n=3
levels <-as.vector(sapply(strsplit(as.character(nmds.df$Taxonomy),"[;] "), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))

nmds.df$phylum <-as.factor(sapply(levels, function(x) paste(x[2:2], collapse=";")))
nmds.df$class <-as.factor(sapply(levels, function(x) paste(x[3:3], collapse=";"))) 
nmds.df$order <-as.factor(sapply(levels, function(x) paste(x[4:4], collapse=";")))
nmds.df$family <-as.factor(sapply(levels, function(x) paste(x[5:5], collapse=";")))
nmds.df$genus <-as.factor(sapply(levels, function(x) paste(x[6:6], collapse=";")))
nmds.df$species <-as.factor(sapply(levels, function(x) paste(x[7:7], collapse=";")))
#
#groupvar <- as.character(groupvar)
#hulls.df = do.call(rbind, lapply(unique(groupvar), function(i) {
#data = dis_nmds$points[which(groupvar == i), ]
#data.frame(groupvar=i,data[chull(data[, 1], data[, 2]),],row.names = NULL)}))

#df <- subset(nmds.df, class=="c__Clostridia")
pdf("NMDS_hq_mash_class+family.pdf",height=15,width=15)
ggplot(nmds.df) +
geom_point(aes(NMDS1,NMDS2, color=factor(class)))+
  geom_text_repel(aes(NMDS1,NMDS2,label = family),size=2)
dev.off()


#species instead
# Read in the distance file, aaf.dist for example
distance <-as.dist(df.data.species)
# Turn it into a matrix
colnames(distance) <- rownames(distance)
dis_matrix <- as.matrix(distance)
# Collaps into two dimension using NMDS argorithm
dis_nmds <- isoMDS(dis_matrix)
# Plot 
points <- as.data.frame(dis_nmds$points)
colnames(points) <- c('NMDS1','NMDS2')
#combine with metadata
nmds.df <- merge(points, metadata, by.x="row.names",by.y="Bin")
#get taxa names
levels <-as.vector(sapply(strsplit(as.character(nmds.df$Taxonomy),"[;] "), function(x) list(trimws(x))))
nmds.df$phylum <-as.factor(sapply(levels, function(x) paste(x[2:2], collapse=";")))
nmds.df$class <-as.factor(sapply(levels, function(x) paste(x[3:3], collapse=";")))
nmds.df$order <-as.factor(sapply(levels, function(x) paste(x[4:4], collapse=";")))
nmds.df$family <-as.factor(sapply(levels, function(x) paste(x[5:5], collapse=";")))
nmds.df$genus <-as.factor(sapply(levels, function(x) paste(x[6:6], collapse=";")))
nmds.df$species <-as.factor(sapply(levels, function(x) paste(x[7:7], collapse=";")))
pdf("NMDS_hq_species_mash_class+family.pdf",height=15,width=15)
ggplot(nmds.df) +
geom_point(aes(NMDS1,NMDS2, color=factor(class)))+
  geom_text_repel(aes(NMDS1,NMDS2,label = species),size=2)
dev.off()


################################booostrapped
#library(pvclust)
#library(vegan)

#result <- pvclust(df.data.krakencont, method.dist="cor", method.hclust="average", nboot=1000)
#pvrect(result, alpha=0.95)

#pdf("/workdir/users/agk85/CDC2/bins/hq_bins_krakencont_heatmap_mash.pdf", height=100, width=100)
#heatmap.2(df.data.krakencont, labRow=labels,labCol=labels,
#trace="none", cexRow=.8,cexCol=0.8,lhei=c(.5, 10),lwid=c(.5,10),margins=c(60,60))
#dev.off()
#
#
#require("RColorBrewer")
#myCol <- colorRampPalette(c("dodgerblue", "black", "yellow"))(100)
#myBreaks <- seq(-2, 2, length.out=101)

#labels = paste(df.metadata.krakencont$Taxonomy, " ", df.metadata.krakencont$Bin,sep="")
#pdf("/workdir/users/agk85/CDC2/bins/hq_bins_krakencont_heatmap_mash2.pdf", height=100, width=100)
#heatmap.2(df.data.krakencont,
#  Rowv=as.dendrogram(result$hclust),
#  Colv=as.dendrogram(result$hclust),
#  #col=myCol,
#  breaks=myBreaks, main="Title",  key=T,  keysize=1.0,scale="none",
#  density.info="none",  reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
#  trace="none",  cexRow=0.2,  cexCol=0.8)
#  #distfun=function(x) dist(x, method="euclidean"),
#  #hclustfun=function(x) hclust(x, method="ward.D2"))
#dev.off()













