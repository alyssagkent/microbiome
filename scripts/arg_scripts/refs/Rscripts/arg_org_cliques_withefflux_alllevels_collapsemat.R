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
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n<-14
mech.col<-col_vector[59:(n+59)]
pie(rep(1,n), col=mech.col)
n <- 7
kewl<-col_vector[1:(n+1)]
pie(rep(1,n), col=kewl)
n<-7
phylum.col<-col_vector[47:(n+47)]
pie(rep(1,n), col=phylum.col)
###################################################################
setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques")
depth = 1
#import the cluster info
arginfo= read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_v_samp_95_95_names_mech.txt", sep="\t",header = T)
indata = read.csv(paste("arg_org_hic_cliques_95_98_",depth, "_2.tbl",sep=""), header=T, sep="\t",row.names=1)
merged =merge(arginfo, indata, by="Cluster", all.y=T)
sub2<- merged[order(merged$ARG_name,merged$Sub_mechanism, merged$Cluster),] #indata$Sample, do first if want to reorder

#remove the efflux genes
sub3 <- sub2
#remove all of the other columns
sub3$Protein<-NULL
sub3$Name<-NULL
sub3$CARD<-NULL
sub3$Resfams<-NULL
sub3$Mechanism<-NULL
sub3$Sub_mechanism<-NULL
sub4 <- sub3

connectdf = sub4[5:ncol(sub4)]
connectsum <- rowSums(connectdf)
#everything will either be normal or multitaxa
sub5 <- subset(sub4,connectsum>0)

#row colors
patient.palette <- colorRampPalette(c("red",'blue',"orange","gray","purple","forestgreen","navy"))
sub5$patients <- as.vector(sapply(strsplit(as.character(sub5$Sample),"-"), function(x) x[[1]]))
patient.col=patient.palette(length(unique(factor(sub5$patients))))[factor(sub5$patients)]
######what about if you aggregate sub5 down to patients:
a = subset(arginfo, select=c("Cluster","Sub_mechanism"))
sub5_a <- merge(sub5, a, by="Cluster")
sub6 <- sub5_a[order(sub5_a$Sub_mechanism),]
sub6$Sample<-NULL
#sub6$Sub_mechanism<-NULL
sub6$Top_ARG<-NULL
sub6$ARG_name<-as.character(sub6$ARG_name)

############################takes awhile for d=2
df_melt <- melt(sub6, id = c("Sub_mechanism","patients", "Cluster", "ARG_name"))
df_melted <- aggregate(value ~ ., max, data=df_melt)

for (lev in 2:7){
levels <-as.vector(sapply(strsplit(as.character(df_melted$variable),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[lev:lev], collapse=";")))
df_melted$fulllevel <- as.factor(sapply(levels, function(x) paste(x[1:lev], collapse=";")))
#aggregate at X level and 
sub7_pre <-aggregate(value ~ patients +  Cluster+ ARG_name + fulllevel, data=df_melted, FUN=sum)
sub7_pre$value <- (sub7_pre$value>0)+0
sub7 <- dcast(sub7_pre, patients+Cluster +ARG_name ~ fulllevel, sum)
connectdf <-sub7[,4:ncol(sub7)]

levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf)),";"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
connectdf_t <- data.frame(t(connectdf))
connectdf_t$groupvar <- groupvar
connectdf_speciesonly <- connectdf_t  #subset(connectdf_t, groupvar!="s__.")
connectdf_speciesonly$groupvar <- NULL
r = rowSums(connectdf_speciesonly)
c = colSums(connectdf_speciesonly)
cdf_speciesonly <-connectdf_speciesonly[,c>0]
sub7_metadata <-sub7[c>0,1:3]
sub7_speciesonly <- cbind(sub7_metadata, t(cdf_speciesonly))

pat_mats<-list()
pats <- unique(sub7_speciesonly$patients)
for (i in seq(length(pats))){
	print(i)
	patient <- pats[i]
	print(patient)
	sub9 <- subset(sub7_speciesonly, patients==patient)
	sub9$ARG_name <- NULL
	sub9$patients <- NULL
	rownames(sub9)<-sub9$Cluster
	sub9$Cluster<-NULL
	sub10<-as.data.frame(t(sub9[,0:(ncol(sub9))]))
	sub11 <- sub10[rowSums(sub10)>0,colSums(sub10)>1]
	x <-sub11
	x <- apply(x, 2,  function(x) as.numeric(x > 0))  #recode as 0/1
	v <- x %*% t(x)                                   #the magic matrix 
	diag(v) <- 0                                      #repalce diagonal
	dimnames(v) <- list(rownames(sub11), rownames(sub11))                #name the dimensions
	v[is.na(v)]<-0
	pat_mats[[i]]<-v
}
allargs_pat_mats <-pat_mats
##############################
m1 <- melt(pat_mats[[1]])
m2 <- melt(pat_mats[[2]])
m3 <- melt(pat_mats[[3]])
m4 <- melt(pat_mats[[4]])
m5 <- melt(pat_mats[[5]])
m6 <- melt(pat_mats[[6]])
m7 <- melt(pat_mats[[7]])

longdf <- rbind(m1,m2,m3,m4,m5,m6,m7)
longdf.agg <- aggregate(value ~., data =longdf,sum)
longdf.link <- subset(longdf.agg, value>0)
el=as.matrix(longdf.link)
g=graph.edgelist(el[,1:2])
E(g)$weight=as.numeric(el[,3])
g_simp1<-as.undirected(g, "collapse")
g_simp <- simplify( g_simp1, remove.multiple=T, edge.attr.comb=c(weight="sum") )
E(g_simp)$weight<-E(g_simp)$weight/2

levels <-as.vector(sapply(strsplit(as.character(V(g_simp)$name),";"), function(x) list(trimws(x))))
V(g_simp)$groupvar <- as.character(as.factor(sapply(levels, function(x) paste(x[5:5], collapse=";"))))
edgelist <-E(g_simp)[inc(V(g_simp)["f__Enterobacteriaceae" ==groupvar])]
E(g_simp)$enterocol <- ifelse(E(g_simp) %in% edgelist, "red","black")

pat_graphs = list()
for (i in seq(length(pats))){
pat.link <- as.matrix(subset(melt(pat_mats[[i]]), value>0))
g.pat <- graph.edgelist(pat.link[,1:2])
E(g.pat)$weight=as.numeric(as.matrix(pat.link)[,3])
g.pat.sim<-as.undirected(g.pat, "collapse")
g.pat.simp <- simplify(g.pat.sim, remove.multiple=T, edge.attr.comb=c(weight="sum") )
E(g.pat.simp)$weight<-E(g.pat.simp)$weight/2
pat_graphs[[i]] = g.pat.simp
}

#create an empty graph without edges
g_simp_noedges<-delete_edges(g_simp, E(g_simp))
g_simp_species <- g_simp #delete_vertices(g_simp, V(g_simp)[species=='s__.'])
g_simp_species_noedges<-delete_edges(g_simp_species, E(g_simp_species))

n=2
levels <-as.vector(sapply(strsplit(as.character(V(g_simp_species_noedges)$name),";"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
familydown <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
familycolor<-col_vector[1:(length(unique(factor(groupvar))))][factor(groupvar)]

n=2
levels <-as.vector(sapply(strsplit(as.character(V(g_simp_species_noedges)$name),";"), function(x) list(trimws(x))))
phyname <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
phylumcolor<-arg.palette(length(unique(factor(phyname))))[factor(phyname)]

V(g_simp_species_noedges)$phylumcolor<- phylumcolor
V(g_simp_species_noedges)$phyname<- as.character(phyname)
V(g_simp_species_noedges)$familydown <- as.character(familydown)
V(g_simp_species_noedges)$familycolor<- familycolor
V(g_simp_species_noedges)$entero <- ifelse(groupvar== "f__Enterobacteriaceae","Enterobacteriaceae","Not Enterobacteriaceae")
V(g_simp_species_noedges)$enterocol <- ifelse(groupvar== "f__Enterobacteriaceae","orange","blue")
graph.to.plot<-g_simp_species_noedges

#for legend
phynamecolor <- data.frame(cbind(as.character(phyname), phylumcolor,V(g_simp_species)$name))
colnames(phynamecolor )<-c("Phylum","Phylumcolor","Fullname")
ppf <- phynamecolor [order(phynamecolor$Fullname),]
ppf$Fullname<-NULL
pp<- unique(ppf)
pp$Phylum<- as.character(pp$Phylum)
pp$Phylumcolor<- as.character(pp$Phylumcolor)
famnamecolor <- data.frame(cbind(as.character(groupvar), familycolor,V(g_simp_species)$name))
colnames(famnamecolor)<-c("Family","Familycolor","Fullname")
fff <- famnamecolor[order(famnamecolor$Fullname),]
fff$Fullname<-NULL
ff<- unique(fff)
ff$Family<- as.character(ff$Family)
ff$Familycolor<- as.character(ff$Familycolor)

#get each patients graph merged with the nodes of everything
g1<-graph.to.plot+pat_graphs[[1]]
g2<-graph.to.plot+pat_graphs[[2]]
g3<-graph.to.plot+pat_graphs[[3]]
g4<-graph.to.plot+pat_graphs[[4]]
g5<-graph.to.plot+pat_graphs[[5]]
g6<-graph.to.plot+pat_graphs[[6]]
g7<-graph.to.plot+pat_graphs[[7]]
patient_graphs=list(g1,g2,g3,g4,g5,g6,g7)

######################plotting the total group as a heatmap
combo_orgorg<-as.matrix(as_adjacency_matrix(g_simp_species,type="both",names=T,attr="weight"))
combo_orgorg<-combo_orgorg[ order(rownames(combo_orgorg)), order(colnames(combo_orgorg))]
n=lev-1
levels <-as.vector(sapply(strsplit(as.character(colnames(combo_orgorg)),";"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:lev], collapse=";")))
colors <- colorRampPalette(brewer.pal(9,"Set1"))
taxa.col = colors(length(unique(groupvar)))
csc=cbind(taxa.col[groupvar])
a <- unique(cbind(as.character(groupvar), csc))
map.palette  <- c("white", colorRampPalette(c("blue", "red"))(n=max(combo_orgorg)))

#make each patient it's individual heatmap
##############plot them all together
matoi <- as.matrix(as_adjacency_matrix(g1,type="both",names=T,attr="weight_2"))
g1_mat <-matoi[order(rownames(matoi)),order(colnames(matoi))]
matoi <- as.matrix(as_adjacency_matrix(g2,type="both",names=T,attr="weight_2"))
g2_mat <-matoi[order(rownames(matoi)),order(colnames(matoi))]
matoi <- as.matrix(as_adjacency_matrix(g3,type="both",names=T,attr="weight_2"))
g3_mat <-matoi[order(rownames(matoi)),order(colnames(matoi))]
matoi <- as.matrix(as_adjacency_matrix(g4,type="both",names=T,attr="weight_2"))
g4_mat <-matoi[order(rownames(matoi)),order(colnames(matoi))]
matoi <- as.matrix(as_adjacency_matrix(g5,type="both",names=T,attr="weight_2"))
g5_mat <-matoi[order(rownames(matoi)),order(colnames(matoi))]
matoi <- as.matrix(as_adjacency_matrix(g6,type="both",names=T,attr="weight_2"))
g6_mat <-matoi[order(rownames(matoi)),order(colnames(matoi))]
matoi <- as.matrix(as_adjacency_matrix(g7,type="both",names=T,attr="weight_2"))
g7_mat <-matoi[order(rownames(matoi)),order(colnames(matoi))]
#####################################################################################################
mybreaks <- seq(0, max(combo_orgorg), length.out=max(combo_orgorg)+2)

pdf(paste("arg_heatmap_org_org_none_separated_", lev, "_withefflux_alllevels_depth1.pdf",sep=""),height=30,width=30)
par(mfrow=c(8,1))

heatmap.3(combo_orgorg,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,breaks=mybreaks,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

heatmap.3(g1_mat,na.rm = TRUE, scale="none", hclustfun=myclust, 
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none", breaks=mybreaks,
main="B314",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

heatmap.3(g2_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",breaks=mybreaks,
main="B316",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

heatmap.3(g3_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",breaks=mybreaks,
main="B320",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

heatmap.3(g4_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",breaks=mybreaks,
main="B331",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

heatmap.3(g5_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",breaks=mybreaks,
main="B335",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

heatmap.3(g6_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",breaks=mybreaks,
main="B357",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

heatmap.3(g7_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",breaks=mybreaks,
main="B370",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=1)

dev.off()
}
