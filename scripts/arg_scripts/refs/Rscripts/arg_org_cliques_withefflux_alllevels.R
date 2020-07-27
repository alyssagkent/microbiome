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
df_melt <- melt(sub6, id = c("patients", "Cluster", "ARG_name","Sub_mechanism"))
#cast with the patients
sub7 <- dcast(df_melt, Sub_mechanism + Cluster+ patients + ARG_name ~ variable, sum)

connectdf <-sub7[,5:ncol(sub7)]
connectdf <- connectdf[,colSums(connectdf)>0]
connectdf <-(connectdf>0) + 0
sub7_metadata <- sub7[,1:4]
sub7_presence <-cbind(sub7_metadata, connectdf)
sub7_melt <- melt(sub7_presence, id = c("patients", "Cluster", "ARG_name", "Sub_mechanism"))
#cast with the patients
sub7_aggregate <- dcast(sub7_melt, Sub_mechanism + Cluster + ARG_name ~ variable, sum)

dim(sub7_aggregate)
connectdf <-sub7_aggregate[,4:ncol(sub7_aggregate)]
n=7
levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
connectdf_t <- data.frame(t(connectdf))
connectdf_t$groupvar <- groupvar
connectdf_speciesonly <- connectdf_t #subset(connectdf_t, groupvar!="s__.")
connectdf_speciesonly$groupvar <- NULL
r = rowSums(connectdf_speciesonly)
c = colSums(connectdf_speciesonly)
cdf_speciesonly <-connectdf_speciesonly[r>0,c>0]
sub7_speciesonly <-sub7_aggregate[c>0,] 
subdata = data.matrix(cdf_speciesonly)
n=5
levels <-as.vector(sapply(strsplit(as.character(rownames(subdata)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))

colors <- c("white","lightblue","blue","purple","yellow","orange","red","pink")
heat.col = colors[0:max(subdata)+1]
n <- max(subdata)
namelist = c("not connected","connected")
names <- namelist[0:max(subdata)+1]
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}
taxa.colors <- colorRampPalette(brewer.pal(9,"Set1"))
taxa.col = taxa.colors(length(unique(groupvar)))
csc=cbind(taxa.col[groupvar])

connectdf <-sub7[,5:ncol(sub7)]
cdf <- connectdf[rowSums(connectdf)>1,colSums(connectdf)>0]
sub8_metadata <- sub7[rowSums(connectdf)>1,1:4]
connectdf <-(cdf>0) + 0
sub8_presence <-cbind(sub8_metadata, connectdf)
sub8_melt <- melt(sub8_presence, id = c("Sub_mechanism","patients", "Cluster", "ARG_name"))
#cast to aggregate on the patients
sub8_aggregate <- dcast(sub8_melt, Cluster + ARG_name ~ variable, sum)

dim(sub8_aggregate)

connectdf <-sub8_aggregate[,3:ncol(sub8_aggregate)]
n=7
levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
connectdf_t <- data.frame(t(connectdf))
connectdf_t$groupvar <- groupvar
connectdf_speciesonly <- connectdf_t  #subset(connectdf_t, groupvar!="s__.")
connectdf_speciesonly$groupvar <- NULL
r = rowSums(connectdf_speciesonly)
c = colSums(connectdf_speciesonly)
cdf_speciesonly <-connectdf_speciesonly[r>0,c>0]
sub8_speciesonly <-sub8_aggregate[c>0,] 
subdata = data.matrix(cdf_speciesonly)
n=5
levels <-as.vector(sapply(strsplit(as.character(rownames(subdata)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))

colors <- c("white","purple","blue","lightblue","yellow","orange","red","pink")
heat.col = colors[0:max(subdata)+1]
n <- max(subdata)
namelist = c("not connected","connected")
names <- namelist[0:max(subdata)+1]
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}
taxa.colors <- colorRampPalette(brewer.pal(9,"Set1"))
taxa.col = taxa.colors(length(unique(as.character(groupvar))))
csc=cbind(taxa.col[groupvar])
d<-data.frame(unique(cbind(as.character(groupvar),csc)))
colnames(d)<- c("Taxa","Taxacolor")

#row colors
#fixing the mechanism colors
merged <- merge(sub8_speciesonly, arginfo, by="Cluster")
mechs = unique(arginfo$Sub_mechanism)[order(unique(arginfo$Sub_mechanism))]
mcol = mech.col[1:length(levels(merged$Sub_mechanism))] 
cols = setNames(mcol, levels(merged$Sub_mechanism))
k = cols[merged$Sub_mechanism]
kmat = t(as.matrix(k))
kuniq = unique(k)
kuniqnames = unique(names(k))

pdf(paste("heatmap_arg_org_multitaxa_depth1_withefflux_alllevels_none.pdf",sep=""), height=60, width=60,useDingbats=FALSE)
heatmap.3(t(subdata), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="none",	Rowv=F,Colv=F,key=F,symbreaks=FALSE, margins=c(30,30),
	symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
	labRow=paste("Cluster",sub8_aggregate$Cluster,sub8_aggregate$ARG_name,sep="_"),	cexCol=.8,
	col=heat.col, cexRow=.4,RowSideColors=kmat,ColSideColors=csc)
legend("topleft",cex=5,title="Patients",legend=c("None","1","2","3","4","5","6"),fill=colors)
legend("left",cex=5,title="Gene Resistance Mechanisms",legend=kuniqnames,fill=kuniq)
legend("topright",cex=2,title="Family",legend=as.character(d$Taxa),fill=as.character(d$Taxacolor))
dev.off()

pdf(paste("heatmap_arg_org_multitaxa_depth1_withefflux_alllevels_row.pdf",sep=""), height=60, width=60,useDingbats=FALSE)
heatmap.3(t(subdata), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="row",	Rowv=T,Colv=F,key=F,symbreaks=FALSE, margins=c(30,30),
	symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
	labRow=paste("Cluster",sub8_aggregate$Cluster,sub8_aggregate$ARG_name,sep="_"),	cexCol=.8,
	col=heat.col, cexRow=.4,RowSideColors=kmat,ColSideColors=csc)
legend("topleft",cex=5,title="Patients",legend=c("None","1","2","3","4","5","6"),fill=colors)
legend("left",cex=5,title="Gene Resistance Mechanisms",legend=kuniqnames,fill=kuniq)
legend("topright",cex=2,title="Family",legend=as.character(d$Taxa),fill=as.character(d$Taxacolor))
dev.off()

pdf(paste("heatmap_arg_org_multitaxa_depth1_withefflux_alllevels_both.pdf",sep=""), height=60, width=60,useDingbats=FALSE)
heatmap.3(t(subdata), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="both",	Rowv=T,Colv=T,key=F,symbreaks=FALSE, margins=c(30,30),
	symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
	labRow=paste("Cluster",sub8_aggregate$Cluster,sub8_aggregate$ARG_name,sep="_"),	cexCol=.8,
	col=heat.col, cexRow=.4,RowSideColors=kmat,ColSideColors=csc)
legend("topleft",cex=5,title="Patients",legend=c("None","1","2","3","4","5","6"),fill=colors)
legend("left",cex=5,title="Gene Resistance Mechanisms",legend=kuniqnames,fill=kuniq)
legend("topright",cex=2,title="Family",legend=as.character(d$Taxa),fill=as.character(d$Taxacolor))
dev.off()
#############################################################
#multipatient
df_melt <- melt(sub6, id = c("Sub_mechanism","patients", "ARG_name","Cluster"))
colnames(df_melt)<- c("Sub_mechanism","Patient", "ARG_name", "Cluster",  "variable", "value" )
sub7 <- dcast(df_melt, Sub_mechanism +Patient+Cluster + ARG_name~ variable, sum)
sub7_metadata <- sub7[1:4]
connectdf <-sub7[,5:ncol(sub7)]
connectdf <- connectdf[,colSums(connectdf)>0]
#make binary
connectdf <-(connectdf>0) + 0
#this tells you which rows have at least one taxa
sub7_metadata$Patcount <- as.numeric(rowSums(connectdf)>0)
sub7_presence <-cbind(sub7_metadata, connectdf)
sub7_melt <- melt(sub7_presence, id = c("Sub_mechanism","Patient", "Cluster","ARG_name","Patcount"))
sub7_aggregate<- dcast(sub7_melt, Sub_mechanism + Cluster +  ARG_name~ variable, sum)
sub7_multipat<- dcast(sub7_metadata, Sub_mechanism + Cluster + ARG_name ~ Patcount, sum)
colnames(sub7_multipat)<- c("Sub_mechanism","Cluster","ARG_name","Patcount")

for (i in 0:5){
	print(i)
sub7_multipatient <- cbind(sub7_multipat[sub7_multipat$Patcount>i,], sub7_aggregate[sub7_multipat$Patcount>i,4:ncol(sub7_aggregate)])
connectdf_multipat <-sub7_multipatient[,5:ncol(sub7_multipatient)]
#trim down to species
n=7
levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf_multipat)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
connectdf_multipat_t <- data.frame(t(connectdf_multipat))
connectdf_multipat_t$groupvar <- groupvar
connectdf_multipat_speciesonly <- connectdf_multipat_t #subset(connectdf_multipat_t, groupvar!="s__.")
connectdf_multipat_speciesonly$groupvar <- NULL
r = rowSums(connectdf_multipat_speciesonly)
c = colSums(connectdf_multipat_speciesonly)
cdf_multipat_speciesonly <-connectdf_multipat_speciesonly[r>0,c>0]
sub7_multipat_speciesonly <-sub7_multipatient[c>0,1:4] 
subdata_multipat = data.matrix(cdf_multipat_speciesonly)
n=5
levels <-as.vector(sapply(strsplit(as.character(rownames(subdata_multipat)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))

colors <- c("white","purple","blue","lightblue","yellow","orange","red","pink")
heat.col = colors[0:max(subdata_multipat)+1]
n <- max(subdata_multipat)
namelist = c("not connected","connected")
names <- namelist[0:max(subdata_multipat)+1]
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}
taxa.colors <- colorRampPalette(brewer.pal(9,"Set1"))
taxa.col = taxa.colors(length(unique(groupvar)))
csc=cbind(taxa.col[groupvar])
d<-data.frame(unique(cbind(as.character(groupvar),csc)))
colnames(d)<- c("Taxa","Taxacolor")

#fixing the mechanism colors
mechs = unique(arginfo$Sub_mechanism)[order(unique(arginfo$Sub_mechanism))]
mcol = mech.col[1:length(levels(arginfo$Sub_mechanism))] 
cols = setNames(mcol, levels(arginfo$Sub_mechanism))
k = cols[sub7_multipat_speciesonly$Sub_mechanism]
kmat = t(as.matrix(k))
kuniq = unique(k)
kuniqnames = unique(names(k))

pdf(paste("heatmap_arg_org_",i+1,"+patients_depth1_withefflux_alllevels_none.pdf",sep=""), height=60, width=40,useDingbats=FALSE)
heatmap.3(t(subdata_multipat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="none",	Rowv=F,Colv=F,key=F,symbreaks=FALSE, margins=c(30,30),
	symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
	labRow=paste("Cluster",sub7_multipat_speciesonly$Cluster,sub7_multipat_speciesonly$ARG_name,sep="_"),	cexCol=.8,
	col=heat.col, cexRow=(i+1)/10,RowSideColors=kmat,ColSideColors=csc)
legend("topleft",cex=5,title="Patients",legend=c("None","1","2","3","4","5","6","7"),fill=colors)
legend("left",cex=3,title="Gene Resistance Mechanisms",legend=kuniqnames,fill=kuniq)
legend("topright",cex=1.75,title="Family",legend=as.character(d$Taxa),fill=as.character(d$Taxacolor))
dev.off()
pdf(paste("heatmap_arg_org_",i+1,"+patients_depth1_withefflux_alllevels_row.pdf",sep=""), height=60, width=40,useDingbats=FALSE)
heatmap.3(t(subdata_multipat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="row",	Rowv=T,Colv=F,key=F,symbreaks=FALSE, margins=c(30,30),
	symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
	labRow=paste("Cluster",sub7_multipat_speciesonly$Cluster,sub7_multipat_speciesonly$ARG_name,sep="_"),	cexCol=.8,
	col=heat.col, cexRow=(i+1)/10,RowSideColors=kmat,ColSideColors=csc)
legend("topleft",cex=5,title="Patients",legend=c("None","1","2","3","4","5","6","7"),fill=colors)
legend("left",cex=3,title="Gene Resistance Mechanisms",legend=kuniqnames,fill=kuniq)
legend("topright",cex=1.75,title="Family",legend=as.character(d$Taxa),fill=as.character(d$Taxacolor))
dev.off()
pdf(paste("heatmap_arg_org_",i+1,"+patients_depth1_withefflux_alllevels_both.pdf",sep=""), height=60, width=40,useDingbats=FALSE)
heatmap.3(t(subdata_multipat), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
	dendrogram="both",	Rowv=T,Colv=T,key=F,symbreaks=FALSE, margins=c(30,30),
	symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
	labRow=paste("Cluster",sub7_multipat_speciesonly$Cluster,sub7_multipat_speciesonly$ARG_name,sep="_"),	cexCol=.8,
	col=heat.col, cexRow=(i+1)/10,RowSideColors=kmat,ColSideColors=csc)
legend("topleft",cex=5,title="Patients",legend=c("None","1","2","3","4","5","6","7"),fill=colors)
legend("left",cex=3,title="Gene Resistance Mechanisms",legend=kuniqnames,fill=kuniq)
legend("topright",cex=1.75,title="Family",legend=as.character(d$Taxa),fill=as.character(d$Taxacolor))
dev.off()
}


#histogram of heatmaps patients vs connections
pdf("Histogram_patients_vs_arg_connections.pdf",height=5,width=5)
i = -1
sub7_multipatient <- cbind(sub7_multipat[sub7_multipat$Patcount>i,], sub7_aggregate[sub7_multipat$Patcount>i,4:ncol(sub7_aggregate)])
connectdf_multipat <-sub7_multipatient[,5:ncol(sub7_multipatient)]
values <- as.vector(as.matrix(connectdf_multipat))
val1 <- values[values>0]
val <- data.frame(val1)
ggplot(val)+
xlab("Patients")+
ylab("Organism-Organism Connections")+
geom_histogram(aes(val1), binwidth=1)+
scale_x_continuous(breaks=seq(0,7,1))
dev.off()

##########################end multipatient#################
############################takes awhile for d=2
df_melt <- melt(sub6, id = c("Sub_mechanism","patients", "Cluster", "ARG_name"))
df_melted <- aggregate(value ~ ., max, data=df_melt)
sub7 <- dcast(df_melted, patients+Cluster + ARG_name ~ variable, sum)
connectdf <-sub7[,4:ncol(sub7)]

n=7
levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf)),"[.][.]"), function(x) list(trimws(x))))
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

n=5
levels <-as.vector(sapply(strsplit(as.character(V(g_simp)$name),"[.][.]"), function(x) list(trimws(x))))
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
#n=7
#levels <- as.vector(sapply(strsplit(as.character(V(g_simp_noedges)$name),"[.][.]"), function(x) list(trimws(x))))
#V(g_simp)$species <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
g_simp_species <- g_simp #delete_vertices(g_simp, V(g_simp)[species=='s__.'])
g_simp_species_noedges<-delete_edges(g_simp_species, E(g_simp_species))

n=5
levels <-as.vector(sapply(strsplit(as.character(V(g_simp_species_noedges)$name),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
familydown <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
familycolor<-col_vector[1:(length(unique(factor(groupvar))))][factor(groupvar)]

n=2
levels <-as.vector(sapply(strsplit(as.character(V(g_simp_species_noedges)$name),"[.][.]"), function(x) list(trimws(x))))
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

l=layout_in_circle(g_simp_species,order =order(V(g_simp_species)$name))
node.size= c(10)
###############taxa labels present in the sample
#load in all the taxa that were ever called for B314---long list
taxalist_subbing=function(infilename) {
taxa = read.table(infilename, header=F,row.names=1,sep = "\t")
taxa = rownames(taxa)
taxa = gsub("; ", "..",taxa)
taxa = gsub(";", ".",taxa)
taxa = gsub("[[]", ".",taxa)
taxa = gsub("[]]", ".",taxa)
taxa = gsub("[']", ".",taxa)
taxa_df <- data.frame(taxa, rep(1,length(taxa)))
return(taxa_df)
}

g1<-pat_graphs[[1]]
g2<-pat_graphs[[2]]
g3<-pat_graphs[[3]]
g4<-pat_graphs[[4]]
g5<-pat_graphs[[5]]
g6<-pat_graphs[[6]]
g7<-pat_graphs[[7]]

all_taxa_df <- taxalist_subbing("all.taxa.txt")
g1_taxa_df <- taxalist_subbing("B314.taxa.txt")
g2_taxa_df <- taxalist_subbing("B316.taxa.txt")
g3_taxa_df <- taxalist_subbing("B320.taxa.txt")
g4_taxa_df <- taxalist_subbing("B331.taxa.txt")
g5_taxa_df <- taxalist_subbing("B335.taxa.txt")
g6_taxa_df <- taxalist_subbing("B357.taxa.txt")
g7_taxa_df <- taxalist_subbing("B370.taxa.txt")
V(g1)$present <- V(g1)$name %in% g1_taxa_df$taxa
V(g2)$present <- V(g2)$name %in% g2_taxa_df$taxa
V(g3)$present <- V(g3)$name %in% g3_taxa_df$taxa
V(g4)$present <- V(g4)$name %in% g4_taxa_df$taxa
V(g5)$present <- V(g5)$name %in% g5_taxa_df$taxa
V(g6)$present <- V(g6)$name %in% g6_taxa_df$taxa
V(g7)$present <- V(g7)$name %in% g7_taxa_df$taxa

a = c(V(g1)$present,V(g2)$present,V(g3)$present,V(g4)$present,V(g5)$present,V(g6)$present,V(g7)$present)
write.csv(a,'Taxa_issues_withefflux_alllevels.csv')

#get each patients graph merged with the nodes of everything
g1<-graph.to.plot+pat_graphs[[1]]
g2<-graph.to.plot+pat_graphs[[2]]
g3<-graph.to.plot+pat_graphs[[3]]
g4<-graph.to.plot+pat_graphs[[4]]
g5<-graph.to.plot+pat_graphs[[5]]
g6<-graph.to.plot+pat_graphs[[6]]
g7<-graph.to.plot+pat_graphs[[7]]

write.graph(g1,file="g1_all_arg_graph_withnodes.gml",format="gml")
write.graph(g2,file="g2_all_arg_graph_withnodes.gml",format="gml")
write.graph(g3,file="g3_all_arg_graph_withnodes.gml",format="gml")
write.graph(g4,file="g4_all_arg_graph_withnodes.gml",format="gml")
write.graph(g5,file="g5_all_arg_graph_withnodes.gml",format="gml")
write.graph(g6,file="g6_all_arg_graph_withnodes.gml",format="gml")
write.graph(g7,file="g7_all_arg_graph_withnodes.gml",format="gml")
write.graph(g_simp, file="g_simp_arg_attributes_withnodes.gml",format="gml")


#load in all the taxa that were ever called for B314---long list
all_taxa_df <- taxalist_subbing("all.taxa.txt")
g1_taxa_df <- taxalist_subbing("B314.taxa.txt")
g2_taxa_df <- taxalist_subbing("B316.taxa.txt")
g3_taxa_df <- taxalist_subbing("B320.taxa.txt")
g4_taxa_df <- taxalist_subbing("B331.taxa.txt")
g5_taxa_df <- taxalist_subbing("B335.taxa.txt")
g6_taxa_df <- taxalist_subbing("B357.taxa.txt")
g7_taxa_df <- taxalist_subbing("B370.taxa.txt")
V(g1)$present <- V(g1)$name %in% g1_taxa_df$taxa 
V(g2)$present <- V(g2)$name %in% g2_taxa_df$taxa 
V(g3)$present <- V(g3)$name %in% g3_taxa_df$taxa 
V(g4)$present <- V(g4)$name %in% g4_taxa_df$taxa 
V(g5)$present <- V(g5)$name %in% g5_taxa_df$taxa 
V(g6)$present <- V(g6)$name %in% g6_taxa_df$taxa 
V(g7)$present <- V(g7)$name %in% g7_taxa_df$taxa 

#add all vertices to figure is it too big??

# write.graph(g1,file="g1_all_allargs_graph.gml",format="gml")
# write.graph(g2,file="g2_all_allargs_graph.gml",format="gml")
# write.graph(g3,file="g3_all_allargs_graph.gml",format="gml")
# write.graph(g4,file="g4_all_allargs_graph.gml",format="gml")
# write.graph(g5,file="g5_all_allargs_graph.gml",format="gml")
# write.graph(g6,file="g6_all_allargs_graph.gml",format="gml")
# write.graph(g7,file="g7_all_allargs_graph.gml",format="gml")
# write.graph(graph.to.plot+g1, file="g_simp_noedges_attributes.gml",format="gml")

#write.graph(g_simp+graph.to.plot,file="g_all_allargs_graph.gml",format="gml")

pdf("Org_org_network_circle_phylum_separated_withefflux_alllevels_depth1.pdf",height=50,width=50)
#1-B314, 2-B316, 3-B320, 4-B331, 5-B335, 6-B357, 7-B370
vsize=2
line_scaling = 5
plot(g1,layout=l,vertex.label=NA,vertex.frame.color=V(g1)$enterocol,vertex.size=vsize,vertex.color=V(g1)$phylumcolor,edge.color=kewl[1],edge.width=E(g1)$weight_2/line_scaling)
par(new=TRUE)
plot(g2,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[2],edge.width=E(g2)$weight_2/line_scaling)
par(new=TRUE)
plot(g3,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[3],edge.width=E(g3)$weight_2/line_scaling)
par(new=TRUE)
plot(g4,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[4],edge.width=E(g4)$weight_2/line_scaling)
par(new=TRUE)
plot(g5,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[5],edge.width=E(g5)$weight_2/line_scaling)
par(new=TRUE)
plot(g6,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[6],edge.width=E(g6)$weight_2/line_scaling)
par(new=TRUE)
plot(g7,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[7],edge.width=E(g7)$weight_2/line_scaling)
legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
#legend("bottomright",cex = 4,pch=21,legend=c(ff$Family),col=c(ff$Familycolor),pt.bg=c(ff$Familycolor))
legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
#make sure that g1 still has the max
maxweight = max(E(g1)$weight_2,E(g2)$weight_2,E(g3)$weight_2,E(g4)$weight_2,E(g5)$weight_2,E(g6)$weight_2,E(g7)$weight_2)
legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,maxweight/line_scaling) ,bty="n")
dev.off() 

patient_graphs=list(g1,g2,g3,g4,g5,g6,g7)


#this will get you the phylum graphs for each patient in circle mode with enterolines
for (i in seq(length(pats))){
	patient_graph = patient_graphs[[i]]
	print(pats[i])
	#highlighting the entero connections
	n=5
	levels <-as.vector(sapply(strsplit(as.character(V(patient_graph)$name),"[.][.]"), function(x) list(trimws(x))))
	V(patient_graph)$groupvar <- as.character(as.factor(sapply(levels, function(x) paste(x[5:5], collapse=";"))))
	edgelist <-E(patient_graph)[inc(V(patient_graph)["f__Enterobacteriaceae" ==groupvar])]
	E(patient_graph)$enterocol <- ifelse(E(patient_graph) %in% edgelist, "red",kewl[i])
	V(patient_graph)$phylumcolor2 <- V(patient_graph )$phylumcolor
	V(patient_graph)$phylumcolor2[!V(patient_graph)$present]<-NA
	pdf(paste("Org_org_network_circle_g",i,"_phylum_withefflux_alllevels_enterolines_depth1.pdf",sep=""),height=50,width=50)
	#1-B314, 2-B316, 3-B320, 4-B331, 5-B335, 6-B357, 7-B370
	vsize=2
	plot(patient_graph ,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=V(patient_graph)$phylumcolor2,edge.color=NA,edge.width=E(patient_graph)$weight_2/line_scaling)
	par(new=TRUE)
	plot(patient_graph,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.width=E(patient_graph)$weight_2/line_scaling,edge.color=E(patient_graph)$enterocol)
	legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
	#legend("bottomright",cex = 4,pch=21,legend=c(ff$Family),col=c(ff$Familycolor),pt.bg=c(ff$Familycolor))
	legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
	legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
	maxweight = max(E(patient_graph)$weight_2)
	legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,maxweight/line_scaling) ,bty="n")
		dev.off()
}

#this will get you the phylum graphs for each patient in circle mode
for (i in seq(length(pats))){
	patient_graph = patient_graphs[[i]]
	print(pats[i])
	#highlighting the entero connections
	n=5
	V(patient_graph)$phylumcolor2 <- V(patient_graph )$phylumcolor
	V(patient_graph)$phylumcolor2[!V(patient_graph)$present]<-NA
	pdf(paste("Org_org_network_circle_g",i,"_phylum_withefflux_alllevels_depth1.pdf",sep=""),height=50,width=50)
	#1-B314, 2-B316, 3-B320, 4-B331, 5-B335, 6-B357, 7-B370
	vsize=2
	plot(patient_graph ,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=V(patient_graph)$phylumcolor2,edge.color=NA,edge.width=E(patient_graph)$weight_2/line_scaling)
	par(new=TRUE)
	plot(patient_graph,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.width=E(patient_graph)$weight_2/line_scaling,edge.color=kewl[i])
	legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
	legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
	legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
	maxweight = max(E(patient_graph)$weight_2)
	legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,maxweight/line_scaling) ,bty="n")
	dev.off()
}



####this is a circle figure with phylum coloring of everything
pdf(paste("Org_org_network_circle_all_phylum_withefflux_alllevels_depth1.pdf",sep=""),height=50,width=50)
vsize=2
plot(graph.to.plot,layout=l,vertex.label=NA,vertex.frame.color=V(graph.to.plot)$enterocol,vertex.size=vsize,vertex.color=V(graph.to.plot)$phylumcolor)
par(new=TRUE)
plot(g_simp_species,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.width=E(g_simp_species)$weight/line_scaling)
legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
maxweight = max(E(g_simp_species)$weight)
legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,max(E(g_simp)$weight/line_scaling)) ,bty="n")
dev.off()
#################################################
####this is a circle figure with phylum coloring of everything with enterolines
pdf(paste("Org_org_network_circle_all_phylum_withefflux_alllevels_enterolines_depth1.pdf",sep=""),height=50,width=50)
vsize=2
plot(graph.to.plot,layout=l,vertex.label=NA,vertex.frame.color=V(graph.to.plot)$enterocol,vertex.size=vsize,vertex.color=V(graph.to.plot)$phylumcolor)
par(new=TRUE)
plot(g_simp_species,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.width=E(g_simp_species)$weight/line_scaling,edge.color=E(g_simp_species)$enterocol)
legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
maxweight = max(E(g_simp_species)$weight)
legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,max(E(g_simp_species)$weight/line_scaling)) ,bty="n")
dev.off()


####CHANGE THE LAYOUT BUT RUN DOWN TO THE BOTTOM OF THIS GROUPING
l=layout_with_mds(g_simp_species)
layoutname = "mds"
vsize=2
####this is a forced diagram figure with phylum coloring of everything
pdf(paste("Org_org_network_",layoutname,"_all_phylum_withefflux_alllevels_enterolines_depth1.pdf",sep=""),height=50,width=50)
plot(graph.to.plot,layout=l,vertex.label=NA,vertex.frame.color=V(graph.to.plot)$enterocol,vertex.size=vsize,vertex.color=V(graph.to.plot)$phylumcolor)
par(new=TRUE)
plot(g_simp_species,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.width=E(g_simp_species)$weight/line_scaling, edge.color=E(g_simp_species)$enterocol)
legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
maxweight = max(E(g_simp_species)$weight)
legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,max(E(g_simp_species)$weight/line_scaling)) ,bty="n")
dev.off()

###do individual ones with a network
for (i in seq(length(pats))){
	patient_graph = patient_graphs[[i]]
	print(pats[i])
	pdf(paste("Org_org_network_", layoutname, "_g",i,"_phylum_withefflux_alllevels_depth1.pdf",sep=""),height=50,width=50)
	plot(graph.to.plot,layout=l,vertex.label=NA,vertex.frame.color=V(graph.to.plot)$enterocol,vertex.size=vsize,vertex.color=V(graph.to.plot)$phylumcolor)
	par(new=TRUE)
	plot(patient_graph,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[i],edge.width=E(patient_graph)$weight_2/line_scaling)
	legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
	#legend("bottomright",cex = 4,pch=21,legend=c(ff$Family),col=c(ff$Familycolor),pt.bg=c(ff$Familycolor))
	legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
	legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
	maxweight = max(E(patient_graph)$weight_2)
	legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,maxweight/line_scaling) ,bty="n")
	dev.off()
}

################
pdf(paste("Org_org_network_",layoutname,"_separated_phylum_withefflux_alllevels_depth1.pdf",sep=""),height=50,width=50)
plot(g1,layout=l,vertex.label=NA,vertex.frame.color=V(g1)$enterocol,vertex.size=vsize,vertex.color=V(g1)$phylumcolor,edge.color=kewl[1],edge.width=E(g1)$weight_2/line_scaling)
par(new=TRUE)
plot(g2,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[2],edge.width=E(g2)$weight_2/line_scaling)
par(new=TRUE)
plot(g3,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[3],edge.width=E(g3)$weight_2/line_scaling)
par(new=TRUE)
plot(g4,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[4],edge.width=E(g4)$weight_2/line_scaling)
par(new=TRUE)
plot(g5,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[5],edge.width=E(g5)$weight_2/line_scaling)
par(new=TRUE)
plot(g6,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[6],edge.width=E(g6)$weight_2/line_scaling)
par(new=TRUE)
plot(g7,layout=l,vertex.label=NA,vertex.frame.color=NA,vertex.size=vsize,vertex.color=NA,edge.color=kewl[7],edge.width=E(g7)$weight_2/line_scaling)
legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
#legend("bottomright",cex = 4,pch=21,legend=c(ff$Family),col=c(ff$Familycolor),pt.bg=c(ff$Familycolor))
legend("topleft",cex = 5,pch=1,legend=c("Enterobacteriaceae","Not Enterobacteriaceae"),col=c("orange","blue"))
legend("topright",cex=5,lty=1,lwd=20,legend=c("B314","B316","B320","B331","B335","B357","B370"),col=c(kewl))
maxweight = max(E(g1)$weight_2,E(g2)$weight_2,E(g3)$weight_2,E(g4)$weight_2,E(g5)$weight_2,E(g6)$weight_2,E(g7)$weight_2)
legend("bottomleft", cex =5, legend=c("min: 1",paste("max: ",maxweight, sep="")), lty=1, lwd=c(1,maxweight/line_scaling) ,bty="n")
dev.off()
###DOWN TO HERE
###################################
############################################################################
######################plotting the total group as a heatmap
combo_orgorg<-as.matrix(as_adjacency_matrix(g_simp_species,type="both",names=T,attr="weight"))
combo_orgorg<-combo_orgorg[ order(rownames(combo_orgorg)), order(colnames(combo_orgorg))]

n=3
levels <-as.vector(sapply(strsplit(as.character(colnames(combo_orgorg)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
colors <- colorRampPalette(brewer.pal(9,"Set1"))
taxa.col = colors(length(unique(groupvar)))
csc=cbind(taxa.col[groupvar])
a <- unique(cbind(as.character(groupvar), csc))
map.palette  <- c("white", colorRampPalette(c("blue", "red"))(n=max(combo_orgorg)))


pdf(paste("heatmap_org_org_cliques_depth1_none_all_withefflux_alllevels.pdf",sep=""), height=50, width=50,useDingbats=FALSE)
heatmap.3(combo_orgorg, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist, 
dendrogram="none",	Rowv=F,Colv=F,key=T,symbreaks=FALSE, margins=c(40,40),
symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
labRow=groupvar_ending,	cexCol=.9, col=map.palette,keysize=0.2,
cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
legend("bottomright",legend=c(a[,1]),
fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=2)
dev.off()


#make each patient it's individual heatmap

#g = g1 g2 g3 g4 g5 g6 g7
for (i in seq(length(pats))){
	patient_graph = patient_graphs[[i]]
	print(pats[i])
	combo_orgorg<-as.matrix(as_adjacency_matrix(patient_graph,type="both",names=T,attr="weight_2"))
	combo_orgorg<-combo_orgorg[ order(rownames(combo_orgorg)), order(colnames(combo_orgorg))]
	
	n=3
	levels <-as.vector(sapply(strsplit(as.character(colnames(combo_orgorg)),"[.][.]"), function(x) list(trimws(x))))
	groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
	groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
	colors <- colorRampPalette(brewer.pal(9,"Set1"))
	taxa.col = colors(length(unique(groupvar)))
	csc=cbind(taxa.col[groupvar])
	a <- unique(cbind(as.character(groupvar), csc))
	map.palette  <- c("white", colorRampPalette(c("blue", "red"))(n=max(combo_orgorg)))
	pdf(paste("heatmap_org_org_cliques_depth1_none_all_withefflux_g", i, ".pdf",sep=""), height=50, width=50,useDingbats=FALSE)
	heatmap.3(combo_orgorg, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist,
	dendrogram="none",      Rowv=F,Colv=F,key=T,symbreaks=FALSE, margins=c(40,40),
	symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
	labRow=groupvar_ending, cexCol=.9, col=map.palette,keysize=0.2,
	cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
	legend("bottomright",legend=c(a[,1]),
	fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=2)
	dev.off()

pdf(paste("heatmap_org_org_cliques_depth1_both_all_withefflux_g", i, ".pdf",sep=""), height=50, width=50,useDingbats=FALSE)
        heatmap.3(combo_orgorg, na.rm = TRUE, scale="none", hclustfun=myclust, distfun=mydist,
        dendrogram="both",Rowv=T,Colv=T,key=T,symbreaks=FALSE, margins=c(40,40),
        symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
        labRow=groupvar_ending, cexCol=.9, col=map.palette,keysize=0.2,
        cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)
        legend("bottomright",legend=c(a[,1]),
        fill=c(as.character(a[,2])), border=FALSE, bty="n", y.intersp = 0.7, cex=2)
        dev.off()
}

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

combo_orgorg<-as.matrix(as_adjacency_matrix(g_simp_species,type="both",names=T,attr="weight"))
combo_orgorg<-combo_orgorg[ order(rownames(combo_orgorg)), order(colnames(combo_orgorg))]

n=3
levels <-as.vector(sapply(strsplit(as.character(colnames(combo_orgorg)),"[.][.]"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
colors <- colorRampPalette(brewer.pal(9,"Set1"))
taxa.col = colors(length(unique(groupvar)))
csc=cbind(taxa.col[groupvar])
a <- unique(cbind(as.character(groupvar), csc))
map.palette  <- c("white", colorRampPalette(c("blue", "red"))(n=max(combo_orgorg)))


pdf(paste("heatmap_org_org_cliques_depth1_none_separated_phylum_withefflux_alllevels_depth1.pdf",sep=""),height=30,width=30)
par(mfrow=c(8,1))

heatmap.3(combo_orgorg,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

heatmap.3(g1_mat,na.rm = TRUE, scale="none", hclustfun=myclust, 
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none", 
main="B314",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9, 
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

heatmap.3(g2_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
main="B316",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

heatmap.3(g3_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
main="B320",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

heatmap.3(g4_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
main="B331",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

heatmap.3(g5_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
main="B335",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

heatmap.3(g6_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
main="B357",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

heatmap.3(g7_mat,na.rm = TRUE, scale="none", hclustfun=myclust,
distfun=mydist,dendrogram="none",Rowv=F,Colv=F,key=T,symbreaks=FALSE,
margins=c(40,40),symkey=F, density.info="none", trace="none",
main="B370",labCol=groupvar_ending,labRow=groupvar_ending, cexCol=.9,
col=map.palette,keysize=0.2,cexRow=.9,RowSideColors=t(csc),ColSideColors=csc)

dev.off()
#################################################################
#doing things with samples too
#plots
#1 every sample together converted to a distance matrix based on their connections
	#sub5_a remove info, turn into distance, plot heatmap
#2
#sub_samp<- sub5
#sub_samp$Top_ARG <-NULL
#sub_samp$ARG_name <-NULL
#rownames(sub_samp)<-sub_samp$Cluster
#colnames(sub_samp)<-

















######################################################
