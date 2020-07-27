library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library("devtools")
library(igraph)
set.seed(42)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
range01 <- function(x)(x-min(x))/diff(range(x))
cRamp <- function(x, palette){
  cols <- colorRamp(palette(100))(range01(x))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}
#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
##################################################
#colors
arg.palette <-colorRampPalette(brewer.pal(12,"Set3"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
n <- 32
kewl<-col_vector[2:(n+2)]
pie(rep(1,n), col=kewl)


#########################################
#load graph
###############
setwd("/workdir/users/agk85/CDC/tables/metagenomes4/")
g = read_graph('/workdir/users/agk85/CDC/tables/metagenomes4/all_metagenome_isolate_99_ncol.txt',directed=FALSE, format = c("ncol"))

df = get.adjacency(g, type=c("both"),attr="weight", names=TRUE, sparse=FALSE)
new_df <- df[ order(row.names(df)),order(colnames(df)) ]
isolates = read.csv('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/isolate_comparisons/isolate_organisms.txt',sep="\t", header=T)
isolate_metadata = read.csv('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/isolate_comparisons/isolate_metadata.txt', sep="\t", header=T)
iso_meta <- isolate_metadata[order(isolate_metadata$Name),]
remove <- c("B335-1Bs","B384-1Bs","L117-3s")
isolate_meta <- iso_meta[!iso_meta$Name %in% remove, ]
#isolates vs genes table
geneorg <- read.csv("/workdir/users/agk85/CDC/tables/metagenomes4/genes_vs_isolates.txt",header=T, row.names=1,sep="\t")
geneorgs <- geneorg[order(row.names(geneorg)),]


gene_counts <- read.csv("/workdir/users/agk85/CDC/tables/metagenomes4/genes_counted.txt",header=F,sep="\t")
colnames(gene_counts)<-c("Sample","Count")
n = length(gene_counts$Count)
x = gene_counts$Count
r = matrix(rep(x,each=n),nrow=n)
c = matrix(rep(x, each=n),ncol=n, byrow=TRUE)
genesum = r + c

dfm <- merge(new_df, isolates, by.x="row.names", by.y="Name")
organisms <- dfm$Species
levels <-as.vector(sapply(strsplit(as.character(rownames(new_df)),"-"), function(x) list(trimws(x))))
patients <- as.factor(sapply(levels, function(x) paste(x[1:1], collapse=";")))
#fixing the mechanism colors
orgs = unique(organism)[order(unique(organisms))]
org.col = kewl[1:length(levels(organisms))] 
cols = setNames(org.col, levels(organisms))
k = cols[organisms]
kmat = t(as.matrix(k))
kuniq = unique(k)
kuniqnames = unique(names(k))
#color coding by patient

pats = unique(patients)[order(unique(patients))]
pat.col = kewl[1:length(levels(patients))]
patcols = setNames(pat.col, levels(patients))
k = patcols[patients]
kmat = t(as.matrix(k))
kuniq = unique(k)
kuniqnames = unique(names(k))


my_palette  <- c("white", colorRampPalette(c("blue", "red"))(n=7000))
myclust=function(c) {hclust(c,method="average")}

pdf("Isolate_metagenome_shared_genes.pdf", height=20,width=20)
heatmap.2(new_df,density.info="none",trace="none",hclustfun=myclust,
col=my_palette, Rowv=T, Colv=T, dendrogram="both",scale="row",
RowSideColors=kmat, ColSideColors=t(kmat))
legend("topright",cex=1,title="Isolate Organisms",legend=kuniqnames,fill=kuniq)
dev.off()


pdf("Isolate_metagenome_shared_genes_nocluster.pdf", height=20,width=20)
heatmap.2(new_df,density.info="none",trace="none",hclustfun=myclust,
col=my_palette, Rowv=F, Colv=F, dendrogram="none",
RowSideColors=kmat, ColSideColors=t(kmat))
legend("topright",cex=1,title="Isolate Organisms",legend=kuniqnames,fill=kuniq)
dev.off()




################################
#ok something like a pca or cca

library(vegan)
dfm <- merge(new_df, isolate_metadata, by.x="row.names", by.y="Name")

#sort by the organism name
df_orgsort <- new_df[order(dfm$Species),order(dfm$Species)] 
organisms2 = dfm$Species[order(dfm$Species)]
#fixing the mechanism colors
orgs2 = unique(organisms2)[order(unique(organisms2))]
org.col2 = kewl[1:length(levels(organisms2))]
cols2 = setNames(org.col2, levels(organisms2))
k2 = cols2[organisms2]
kmat2 = t(as.matrix(k2))
kuniq2 = unique(k2)
kuniqnames2 = unique(names(k2))

pdf("Isolate_shared_genes_orgsort.pdf", height=20,width=20)
heatmap.2(df_orgsort,density.info="none",trace="none",hclustfun=myclust,
col=my_palette, Rowv=F, Colv=F, dendrogram="none",
RowSideColors=kmat2, ColSideColors=t(kmat2))
legend("topright",cex=1,title="Isolate Organisms",legend=kuniqnames2,fill=kuniq2)
dev.off()



###ok what if you want to set all the isolates that have the same species to 0 
df.melt <- melt(df)
m1 <- merge(df.melt, isolates, by.x="Var1", by.y="Name") 
m2 <- merge(m1, isolates, by.x="Var2", by.y="Name")
colnames(m2) <- c("Name1","Name2","Genes","Species1","Species2")

ind <- m2$Species1 == m2$Species2
m2[ind, c("Genes")] <- 0
m3 <- m2
m3$Species1 <- NULL
m3$Species2 <- NULL
m4 <- cast(m3, Name1~Name2)
rownames(m4) <- m4$Name1
m4$Name1 <- NULL
m5<-data.matrix(m4)

new_df <- m5[ order(row.names(m5)),order(colnames(m5)) ]
dfm <- merge(new_df, isolates, by.x="row.names", by.y="Name")
organisms <- dfm$Species
#sort by the organism name
df_orgsort <- new_df[order(dfm$Species),order(dfm$Species)]
organisms2 = dfm$Species[order(dfm$Species)]
#fixing the mechanism colors
orgs2 = unique(organisms2)[order(unique(organisms2))]
org.col2 = kewl[1:length(levels(organisms2))]
cols2 = setNames(org.col2, levels(organisms2))
k2 = cols2[organisms2]
kmat2 = t(as.matrix(k2))
kuniq2 = unique(k2)
kuniqnames2 = unique(names(k2))

pdf("Isolate_shared_genes_orgsort_zeroed.pdf", height=20,width=20)
heatmap.2(df_orgsort,density.info="none",trace="none",hclustfun=myclust,
col=my_palette, Rowv=F, Colv=F, dendrogram="none",
RowSideColors=kmat2, ColSideColors=t(kmat2))
legend("topright",cex=1,title="Isolate Organisms",legend=kuniqnames2,fill=kuniq2)
dev.off()


