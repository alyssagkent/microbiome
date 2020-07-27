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

###########
#import data
###########
setwd("/workdir/users/agk85/CDC2/bins")
metaphlan<- read.csv("/workdir/users/agk85/CDC2/metaphlan/cdc/mgm/CDC_mgm_metaphlan.txt",header=T,sep="\t")

taxnames = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
levs = c(2,3,4,5,6,7)
l = 7
taxalev=taxnames[l]
levels <-as.vector(sapply(strsplit(as.character(metaphlan$ID),"\\|"), function(x) length(x)))
metaphlan$levels <- levels
metaphlan_sub <- subset(metaphlan, levels == l)

#this gets a single name from the ID
taxa <-sapply(strsplit(as.character(metaphlan_sub$ID),"__"), `[`, length(strsplit(as.character(metaphlan_sub$ID),"__")[[l+1]]))
rownames(metaphlan_sub)<-taxa

metaphlan_sub$levels <-NULL
################
#colors
################
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 8
kewl<-col_vector[60:(n+60)]
################################
mypalette  <-colorRampPalette(c("red", "blue"))(n=50)
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}

metaphlan_sub$ID <-NULL
Patient<-sapply(strsplit(colnames(metaphlan_sub),"\\."), function(x) x[1])

patients <- factor(Patient)
pat <- data.frame(patients)
dcol = data.frame(patient = levels(patients),color =kewl)
rsc = as.character(dcol$color[match(pat$patients, dcol$patient)])


df4 <- t(data.matrix(metaphlan_sub))

df5.dist <- vegdist(as.matrix(df4), method="jaccard", binary=T)
df5.dist.mat <- as.matrix(df5.dist)
write.table(df5.dist.mat, "sample_jaccard_distance_metaphlan.txt", sep="\t", quote=FALSE)
df5.dist.mat[df5.dist.mat==0]<-NA

pdf("sample_jaccard_distance_metaphlan.pdf", height=20, width=20,useDingbats=FALSE)
heatmap.2((df5.dist.mat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
        dendrogram="both", Rowv=T,Colv=T,key=TRUE,symbreaks=FALSE,
	density.info="none", trace="none", na.color = "white",
	RowSideColor=rsc, ColSideColors=rsc,
	labCol=colnames(df5.dist.mat), labRow=rownames(df5.dist.mat),
	cexCol=2, col=mypalette,keysize=2,key.par = list(cex=2))
dev.off()
###############################################################################33
df5.dist <- vegdist(as.matrix(df4), method="bray")
df5.dist.mat <- as.matrix(df5.dist)
write.table(df5.dist.mat, "sample_bray_distance_metaphlan.txt", sep="\t", quote=FALSE)
df5.dist.mat[df5.dist.mat==0]<-NA

pdf("sample_bray_distance_metaphlan.pdf", height=20, width=20,useDingbats=FALSE)
heatmap.2((df5.dist.mat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
        dendrogram="both", Rowv=T,Colv=T,key=TRUE,symbreaks=FALSE,
        density.info="none", trace="none", na.color = "white",
        RowSideColor=rsc, ColSideColors=rsc,
        labCol=colnames(df5.dist.mat), labRow=rownames(df5.dist.mat),
        cexCol=2, col=mypalette,keysize=2,key.par = list(cex=2))
dev.off()

#one aggregate by patient

dft <- t(metaphlan_sub)
#one.5 convert less than 1% to 0
dft[dft<1]<-0
#aggregate
dflev <-data.frame(dft)
dflev$pat <- factor(Patient)

aggdata <-aggregate(dft, by=list(Patient),
  FUN=sum)
rownames(aggdata)<-aggdata$Group.1
aggdata$Group.1<-NULL
#two convert to a binary matrix
aggdata[aggdata>0]<-1
#compute jaccard
df4 <- (data.matrix(aggdata))

df5.dist <- vegdist(as.matrix(df4), method="jaccard", binary=T)
df5.dist.mat <- as.matrix(df5.dist)
write.table(df5.dist.mat, "patient_jaccard_distance_metaphlan.txt", sep="\t", quote=FALSE)
df5.dist.mat[df5.dist.mat==0]<-NA

pdf("patient_jaccard_distance_metaphlan.pdf", height=20, width=20,useDingbats=FALSE)
heatmap.2((df5.dist.mat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
        dendrogram="both", Rowv=T,Colv=T,key=TRUE,symbreaks=FALSE,
        density.info="none", trace="none", na.color = "white",
        RowSideColor=kewl, ColSideColors=kewl,
        labCol=colnames(df5.dist.mat), labRow=rownames(df5.dist.mat),
        cexCol=2, col=mypalette,keysize=2,key.par = list(cex=2))
dev.off()



