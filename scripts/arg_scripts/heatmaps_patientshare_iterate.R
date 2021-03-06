library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library("devtools")
#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
set.seed(42)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
  theme(panel.grid = element_blank())}
theme_set(theme_nogrid())

args <- commandArgs(trailingOnly = TRUE)
#genetype = 'arg'
#minreads = 2
#minpatients = 4
#inhandle = 'Together_das_2_argtaxa_patientcount.txt'

genetype = args[1]
minreads = args[2]
minpatients = args[3]
inhandle = args[4]

##################################################
#colors
arg.palette <-colorRampPalette(brewer.pal(12,"Set3"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n<-14
mech.col<-col_vector[59:(n+59)]
n <- 12
kewl<-col_vector[60:(n+60)]
###################################################################
setwd(paste("/workdir/users/agk85/CDC2/bins", sep=""))
#import the cluster info
indata = read.table(inhandle, header=F, sep="\t")
header = c("Number_patients","Minreads","Threshtype","Bintype","Level","Cluster","Genetype","Taxonomy")
colnames(indata) <- header
subdata = subset(indata, Threshtype == 'contacts' & Bintype == 'anybin' & Level == 's__')
minreads =indata[1,2]
subdata$Threshtype = NULL
subdata$Bintype = NULL
subdata$Level = NULL
subdata$Minreads = NULL
subdata$Number_patients <- as.numeric(as.character(subdata$Number_patients))

#cast into wide
subdata.mat <- dcast(subdata, Cluster+Genetype  ~Taxonomy,value.var="Number_patients")
subdata.mat[is.na(subdata.mat)] <- 0
rownames(subdata.mat) <- subdata.mat$Cluster
subdata.mat$Cluster <-NULL
gt = subdata.mat$Genetype
subdata.mat$Genetype <- NULL
subdata.mat[is.na(subdata.mat)] <- 0

sm <- data.matrix(subdata.mat)

#make the boolean matrix of greater than equal to minpatinets
sm_bool<- sm>=minpatients

#only take the rows where the rowsum>0
sm_bool_rs <- rowSums(sm_bool)
sm_genemin <- subset(sm, sm_bool_rs>0)

#only take the genetypes where the row had 1+ connections>minpatients
gt_genemin <- gt[sm_bool_rs>0]

colors <- c("white","azure4","darkmagenta","blue","lightblue","aquamarine3","darkgoldenrod2","chocolate2","darkred","hotpink2")
heat.col = colors[0:max(sm_genemin)+1]
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}

#set groupvar as the family name
n=3
levels <-as.vector(sapply(strsplit(as.character(colnames(sm_genemin)),";"), function(x) list(trimws(x))))
groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))

n=5
groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))

n <- max(sm_genemin)
taxa.colors <- colorRampPalette(brewer.pal(9,"Set1"))
taxa.col = taxa.colors(length(unique(as.character(groupvar))))
csc=cbind(taxa.col[groupvar])
d<-data.frame(unique(cbind(as.character(groupvar),csc)))
colnames(d)<- c("Taxa","Taxacolor")

#row colors
#fixing the mechanism colors
mcol = mech.col[1:length(levels(gt_genemin))] 
cols = setNames(mcol, levels(gt_genemin))
k = cols[gt_genemin]
kmat = t(as.matrix(k))
kuniq = unique(k)
kuniqnames = unique(names(k))

#this is the pdf I want
pdf(paste("heatmap_patientshare_",genetype,"_",minreads,"_reads_",minpatients,"_minpatients_gene.pdf",sep=""), height=30, width=40,useDingbats=FALSE)
heatmap.3(sm_genemin, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none",
        dendrogram="row",Rowv=T,Colv=F,key=F,symbreaks=FALSE, margins=c(20,30),
        symkey=F, density.info="none", trace="none", labCol=groupvar_ending,
        labRow=paste("Cluster",rownames(sm_genemin),gt_genemin,sep="_"),cexCol=.8,
        col=heat.col, cexRow=.5,RowSideColors=kmat,ColSideColors=csc)
legend("topleft",cex=5,title="Number of patients",legend=c("None","1","2","3","4","5","6","7","8","9"),fill=colors)
legend("top",cex=3,title="Genetypes",legend=kuniqnames,fill=kuniq)
legend("topright",cex=1.75,title="Class",legend=as.character(d$Taxa),fill=as.character(d$Taxacolor))
dev.off()
