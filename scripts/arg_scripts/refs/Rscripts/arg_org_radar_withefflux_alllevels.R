library(ggplot2)
library(gplots)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library("devtools")
library(igraph)
library(fmsb)

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

#connectdf <-sub7[,5:ncol(sub7)]
#connectdf <- connectdf[,colSums(connectdf)>0]
#connectdf <-(connectdf>0) + 0
#sub7_metadata <- sub7[,1:4]
#sub7_presence <-cbind(sub7_metadata, connectdf)
#sub7_melt <- melt(sub7_presence, id = c("patients", "Cluster", "ARG_name", "Sub_mechanism"))
##cast with the patients
#sub7_aggregate <- dcast(sub7_melt, Sub_mechanism + Cluster + ARG_name ~ variable, sum)

#dim(sub7_aggregate)
#connectdf <-sub7_aggregate[,4:ncol(sub7_aggregate)]
#n=7
#levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf)),"[.][.]"), function(x) list(trimws(x))))
#groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
#groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))
#connectdf_t <- data.frame(t(connectdf))
#connectdf_t$groupvar <- groupvar
#connectdf_speciesonly <- connectdf_t #subset(connectdf_t, groupvar!="s__.")
#connectdf_speciesonly$groupvar <- NULL
#r = rowSums(connectdf_speciesonly)
#c = colSums(connectdf_speciesonly)
#cdf_speciesonly <-connectdf_speciesonly[r>0,c>0]
#sub7_speciesonly <-sub7_aggregate[c>0,] 
#subdata = data.matrix(cdf_speciesonly)
#n=5
#levels <-as.vector(sapply(strsplit(as.character(rownames(subdata)),"[.][.]"), function(x) list(trimws(x))))
#groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
#groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))

#colors <- c("white","lightblue","blue","purple","yellow","orange","red","pink")
#heat.col = colors[0:max(subdata)+1]
#n <- max(subdata)
#namelist = c("not connected","connected")
#names <- namelist[0:max(subdata)+1]
#mydist=function(c) {dist(c)}
#myclust=function(c) {hclust(c,method="average")}
#taxa.colors <- colorRampPalette(brewer.pal(9,"Set1"))
#taxa.col = taxa.colors(length(unique(groupvar)))
#csc=cbind(taxa.col[groupvar])
#
#connectdf <-sub7[,5:ncol(sub7)]
#cdf <- connectdf[rowSums(connectdf)>1,colSums(connectdf)>0]
#sub8_metadata <- sub7[rowSums(connectdf)>1,1:4]
#connectdf <-(cdf>0) + 0
#sub8_presence <-cbind(sub8_metadata, connectdf)
#sub8_melt <- melt(sub8_presence, id = c("Sub_mechanism","patients", "Cluster", "ARG_name"))
##cast to aggregate on the patients
#sub8_aggregate <- dcast(sub8_melt, Cluster + ARG_name ~ variable, sum)

#dim(sub8_aggregate)

#connectdf <-sub8_aggregate[,3:ncol(sub8_aggregate)]
#n=7
#levels <-as.vector(sapply(strsplit(as.character(colnames(connectdf)),"[.][.]"), function(x) list(trimws(x))))
#groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
#connectdf_t <- data.frame(t(connectdf))
#connectdf_t$groupvar <- groupvar
#connectdf_speciesonly <- connectdf_t  #subset(connectdf_t, groupvar!="s__.")
#connectdf_speciesonly$groupvar <- NULL
#r = rowSums(connectdf_speciesonly)
#c = colSums(connectdf_speciesonly)
#cdf_speciesonly <-connectdf_speciesonly[r>0,c>0]
#sub8_speciesonly <-sub8_aggregate[c>0,] 
#subdata = data.matrix(cdf_speciesonly)
#n=5
#levels <-as.vector(sapply(strsplit(as.character(rownames(subdata)),"[.][.]"), function(x) list(trimws(x))))
#groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
#groupvar_ending <- as.factor(sapply(levels, function(x) paste(x[n:7], collapse=";")))

colors <- c("white","purple","blue","lightblue","yellow","orange","red","pink")
#heat.col = colors[0:max(subdata)+1]
#n <- max(subdata)
#namelist = c("not connected","connected")
#names <- namelist[0:max(subdata)+1]
mydist=function(c) {dist(c)}
myclust=function(c) {hclust(c,method="average")}
#taxa.colors <- colorRampPalette(brewer.pal(9,"Set1"))
#taxa.col = taxa.colors(length(unique(as.character(groupvar))))
#csc=cbind(taxa.col[groupvar])
#d<-data.frame(unique(cbind(as.character(groupvar),csc)))
#colnames(d)<- c("Taxa","Taxacolor")

#row colors
#fixing the mechanism colors
#merged <- merge(sub8_speciesonly, arginfo, by="Cluster")
#mechs = unique(arginfo$Sub_mechanism)[order(unique(arginfo$Sub_mechanism))]
#mcol = mech.col[1:length(levels(merged$Sub_mechanism))] 
#cols = setNames(mcol, levels(merged$Sub_mechanism))
#k = cols[merged$Sub_mechanism]
#kmat = t(as.matrix(k))
#kuniq = unique(k)
#kuniqnames = unique(names(k))

#############################################################
#fixing the mechanism colors
mechs = unique(arginfo$Sub_mechanism)[order(unique(arginfo$Sub_mechanism))]
mcol = mech.col[1:length(levels(arginfo$Sub_mechanism))] 
cols = setNames(mcol, levels(arginfo$Sub_mechanism))
#k = cols[sub7_multipat_speciesonly$Sub_mechanism]
#kmat = t(as.matrix(k))
#kuniq = unique(k)
#kuniqnames = unique(names(k))
##########################end multipatient#################
makeTransparent<-function(someColor, alpha)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

df_melt <- melt(sub6, id = c("Sub_mechanism","patients", "Cluster", "ARG_name"))
df_melted <- aggregate(value ~ ., max, data=df_melt)
df_melted$levels <-as.vector(sapply(strsplit(as.character(df_melted$variable),"[.][.]"), function(x) list(trimws(x))))

#df_melted$Family <- as.factor(sapply(df_melted$levels, function(x) paste(x[1:5], collapse=";")))
#df_melted.agg <- aggregate(value ~ Sub_mechanism+patients + Cluster + ARG_name + Family, sum, data=df_melted)
#df_melted.agg$value <-(df_melted.agg$value>0) + 0


#df_melted$Family <- as.factor(sapply(df_melted$levels, function(x) paste(x[1:7], collapse=";")))
#df_melted.agg <- aggregate(value ~ Sub_mechanism+patients + Cluster + ARG_name + Family, sum, data=df_melted)
#df_melted.agg$value <-(df_melted.agg$value>0) + 0



#do something where you aggregate enteros
df_melted$Family <- as.factor(sapply(df_melted$levels, function(x) paste(x[1:5], collapse=";")))
df_melted_enteros<-subset(df_melted, Family == "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae")
df_melted_enteros.agg <- aggregate(value ~ Sub_mechanism+patients + Cluster + ARG_name + Family, sum, data=df_melted_enteros)

#everyone else
df_melted_notenteros<-subset(df_melted, Family != "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae")
df_melted_notenteros$Family <- as.factor(sapply(df_melted_notenteros$levels, function(x) paste(x[1:7], collapse=";")))
#this probably shouldn't do anythign but remove levels and variable
df_melted_notenteros.agg<-aggregate(value ~ Sub_mechanism+patients + Cluster + ARG_name + Family, sum, data=df_melted_notenteros)

df_entero_notenteros <- rbind(df_melted_enteros.agg, df_melted_notenteros.agg)
df_entero_notenteros$value <-(df_entero_notenteros$value>0) + 0


sub7 <- dcast(df_melted.agg, patients+Cluster+ARG_name~Family, sum)
connectdf <-sub7[,4:ncol(sub7)]
connectdf_t <- data.frame(t(connectdf))
sub7_metadata <-sub7[1:3]
sub7_speciesonly <- cbind(sub7_metadata, t(connectdf_t))


sub7 <- dcast(df_entero_notenteros, patients+Cluster+ARG_name~Family, sum)
connectdf <-sub7[,4:ncol(sub7)]
connectdf_t <- data.frame(t(connectdf))
sub7_metadata <-sub7[1:3]
sub7_speciesonly <- cbind(sub7_metadata, t(connectdf_t))



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
	print(length(unique(sub9$Cluster)))
	print(dim(sub9))
	sub9$Cluster<-NULL
	sub10<-as.data.frame(t(sub9))
	sub11 <- sub10 #[rowSums(sub10)>0,colSums(sub10)>1]
	x <-sub11
	x <- apply(x, 2,  function(x) as.numeric(x > 0))  #recode as 0/1
	v <- x %*% t(x)                                   #the magic matrix 
	diag(v) <- 0                                      #repalce diagonal
	dimnames(v) <- list(rownames(sub11), rownames(sub11))                #name the dimensions
	v[is.na(v)]<-0
	pat_mats[[i]]<-v
}
##############################
psummedlist <- list()
for (i in seq(length(pats))){
p <- data.frame(pat_mats[[i]])
n=5
#n=7
levels <-as.vector(sapply(strsplit(as.character(rownames(p)),";"), function(x) list(trimws(x))))
p$groupvar <- as.factor(sapply(levels, function(x) paste(x[n:n], collapse=";")))
penteros = subset(p, groupvar == "f__Enterobacteriaceae")
#penteros= subset(p, groupvar == "s__Escherichia_coli.")
penteros$groupvar <- NULL
psummed = colSums(penteros)
psummedlist[[i]] = psummed
}

makeTransparent<-function(someColor, alpha)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

enterohits <- cbind(psummedlist[[1]],psummedlist[[2]], psummedlist[[3]], psummedlist[[4]], psummedlist[[5]], psummedlist[[6]], psummedlist[[7]])
colnames(enterohits)<-pats
data <- t(enterohits)
data=rbind(rep(max(data),ncol(data)) , rep(0,ncol(data)) , data)
data = data.frame(data)
data = data[ , order(names(data))]
data$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Enterobacteriaceae<-NULL
data$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli.<-NULL
data<- data[,colSums(data[3:9,])>0]


levels <-as.vector(sapply(strsplit(as.character(colnames(data)),"[.]"), function(x) list(trimws(x))))
phyname <- as.factor(sapply(levels, function(x) paste(x[2:2], collapse=";")))
phylumcolor<-arg.palette(length(unique(factor(phyname))))[factor(phyname)]
phynamecolor <- data.frame(cbind(as.character(phyname), phylumcolor,colnames(data)))
colnames(phynamecolor )<-c("Phylum","Phylumcolor","Fullname")
ppf <- phynamecolor [order(phynamecolor$Fullname),]
ppf$Fullname<-NULL
pp<- unique(ppf)
pp$Phylum<- as.character(pp$Phylum)
pp$Phylumcolor<- as.character(pp$Phylumcolor)


#pdf("Giant_radarplot.pdf",height = 50, width=50)
#radarchart(data, pcol=kewl, axistype=1,maxmin=T,plty=1,plwd=5,
#pdensity=30,pfcol=makeTransparent(kewl,30),
#cglcol=phylumcolor, cglty=1, axislabcol="black", caxislabels=seq(0,80,20), cglwd=5)
#legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
#dev.off()

pdf("Giant_radarplot_species_to_entero_nofill.pdf",height = 50, width=50)
radarchart(data, pcol=kewl, axistype=4,maxmin=T,plty=1,plwd=5,vlabels=NULL,
caxislabels = seq(from = min(data), to = max(data), length = 5),
cglcol=phylumcolor, cglty=1, axislabcol="black", cglwd=5)
legend("bottomright",cex = 4,pch=21,legend=c(pp$Phylum),col=c(pp$Phylumcolor),pt.bg=c(pp$Phylumcolor))
dev.off()














