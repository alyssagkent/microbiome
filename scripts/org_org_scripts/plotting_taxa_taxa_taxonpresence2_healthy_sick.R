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
    theme(panel.grid = element_blank())}
theme_set(theme_nogrid())

infile = 'arg_das_2_taxa_taxa_healthysick.txt'
genetype = 'arg'
contacts = 2
line_scaling = 5

args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
contacts = args[2]
genetype = args[3]
line_scaling = as.numeric(args[4])

###########
#import data
###########
setwd("/workdir/users/agk85/CDC2/bins")
df <-read.table(infile,header=T, row.names=1,sep="\t")
################
#colors
################
arg.palette <-colorRampPalette(brewer.pal(12,"Set3"))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 9
kewl<-col_vector[60:(n+60)]
################################


taxalist_subbing=function(infilename,level) {
taxa = read.table(infilename, header=F,row.names=1,sep = "\t")
taxa = rownames(taxa)
taxlev <-as.vector(sapply(strsplit(as.character(taxa),"[;][ ]"), function(x) list(trimws(x))))
taxalev <- as.character(as.factor(sapply(taxlev,function(x) paste(x[level:level], collapse=";"))))
#taxa = gsub("; ", "..",taxalev)
#taxa = gsub(";", ".",taxa)
#taxa = gsub("[[]", ".",taxa)
#taxa = gsub("[]]", ".",taxa)
#taxa = gsub("[']", ".",taxa)
taxa_df <- data.frame(taxalev, rep(1,length(taxa)))
return(taxa_df)
}
df$Genes <- NULL

levels = c('k__','p__','c__')
patients = c('sick','healthy')
for (j in seq(2,length(levels))){
	level = levels[j]
	print(level)
	densities <- vector( length=2)
	for (i in seq(length(patients))){
		patient = patients[i]
		print(patient)
		df_sub <- subset(df, Patient==patient & Level==level)
		c1 = data.frame(as.character(df_sub$Taxa1), as.character(df_sub$Base1))
		colnames(c1)<-c("Taxa","B")
		c2 = data.frame(as.character(df_sub$Taxa2), as.character(df_sub$Base2))
		colnames(c2)<-c("Taxa","B")
		c3 = rbind(c1,c2)
		c4 = unique(c3)
		#get the graph
		a <- subset(df_sub, select=c(Taxa1,Taxa2,Genecount))
		colnames(a) <- c("ref1","ref2","weight")
		g<-graph_from_data_frame(a, directed=FALSE)
		g1 <- simplify( g, remove.multiple=T,remove.loops=F, edge.attr.comb=c(weight="sum"))
		E(g1)$weight <- E(g1)$weight
		d = edge_density(g1)
		densities[i] = d
		if (patient == 'sick'){
                        g1_noedges<-delete_edges(g1, E(g1))
                        n=2
                        taxlev <-as.vector(sapply(strsplit(as.character(V(g1)$name),"[;][ ]"), function(x) list(trimws(x))))
                        V(g1_noedges)$taxname <- as.character(as.factor(sapply(taxlev,function(x) paste(x[j:j], collapse=";"))))
                        l=layout_in_circle(g1,order =order(V(g1)$name))
                        #the all group should have the maximum line weight---otherwise something went wrong!
                        maxweight = max(E(g1)$weight)
                }
		patientmaxweight = max(E(g1)$weight)
		#add the noedges plus the graph of interest together
		gtp <- g1_noedges+g1
		pdf(paste("org_org_figures/Org_org_", genetype, "_",patient,"_",contacts, "__",level,"circle_hs.pdf",sep=""),height=10,width=10)
		vsize=2
		plot(gtp,layout=l,vertex.label.cex=1, vertex.label.dist=0.5, vertex.color='black',vertex.label=V(gtp)$taxname, vertex.size=vsize,edge.width=E(gtp)$weight_2/line_scaling)
		legend("bottomleft", cex =1, legend=c("min: 1",paste("half: ", patientmaxweight/2, sep=""),paste("max: ",patientmaxweight, sep="")), lty=1, lwd=c(1,patientmaxweight/line_scaling/2,patientmaxweight/line_scaling) ,bty="n")
		par(new=TRUE)
		dev.off()
	}
	print(level)
	print(densities)	
}


