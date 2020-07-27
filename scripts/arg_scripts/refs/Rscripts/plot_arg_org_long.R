library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())

setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/")
indata = read.csv("arg_org_gene_hic_95_98_1_2.tbl",header = T,sep="\t", row.names=NULL)
indata$Organism <- NULL
colnames(indata)[1:3] <- c("Cluster", "ARG_name","Organism")
sub = subset(indata, ARG_name == 'AAC3')
sub = subset(sub, 
subdata = data.matrix(sub[,4:ncol(sub)])
colors <- colorRampPalette(c("black","white","red","blue","orange"))
heat.col = colors(max(subdata)+1)
n <- max(mat)
prefix <- 0:n
suffix <- "coded"
names=paste(prefix, suffix, sep=" ")

mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="single")}
pdf(paste("heatmap_arg_org_",patient,"_", pid, "_", mapid, "_", depth, "_", hicweight, "_one.pdf",sep=""), height=10, width=15,useDingbats=FALSE)
heatmap.2(subdata, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="col",	Rowv=F,Colv=T,key=F,symbreaks=FALSE, margins=c(10,15),
	symkey=F, density.info="none", trace="none", labCol=colnames(mat),
	labRow=rownames(mat),	cexRow=.4, cexCol=.4, col=heat.col)
legend("left",legend=c(as.character(d$levels),as.character(a$firstword), names),
	fill=c(as.character(d$rsc),as.character(a$csc),heat.col), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

dev.off()



