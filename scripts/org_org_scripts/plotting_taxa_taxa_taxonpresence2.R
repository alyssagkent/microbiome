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


levels = c('k__','p__','c__','o__','f__','g__', 's__')
patients = c('all','B314','B316','B320','B331','B335','B357','B370','US3','US8')
for (j in seq(2,length(levels))){
	level = levels[j]
	print(level)
	densities <- vector( length=9)
	for (i in seq(length(patients))){
		patient = patients[i]
		print(patient)
		df_sub <- subset(df, Patient==patient & Level==level)
		c1 = data.frame(as.character(df_sub$Taxa1), as.character(df_sub$Color1))
		colnames(c1)<-c("Taxa","Color")
		c2 = data.frame(as.character(df_sub$Taxa2), as.character(df_sub$Color2))
		colnames(c2)<-c("Taxa","Color")
		c3 = rbind(c1,c2)
		c4 = unique(c3)
		#get the graph
		a <- subset(df_sub, select=c(Taxa1,Taxa2,Genes))
		colnames(a) <- c("ref1","ref2","weight")
		g<-graph_from_data_frame(a, directed=FALSE)
		g1 <- simplify( g, remove.multiple=T, remove.loops=F,edge.attr.comb=c(weight="sum"))
		d = edge_density(g1)
		densities[i] = d
		if (patient == 'all'){
			g1_noedges<-delete_edges(g1, E(g1))
			V(g1_noedges)$color <- c4$Color
			n=2
			taxlev <-as.vector(sapply(strsplit(as.character(V(g1)$name),"[;][ ]"), function(x) list(trimws(x))))
			allleg <- c4
			allleg$phyla <- as.character(as.factor(sapply(taxlev, function(x) paste(x[n:n], collapse=";"))))
			allleg$color<-arg.palette(length(unique(factor(allleg$phyla))))[factor(allleg$phyla)]
			V(g1_noedges)$taxname <- as.character(as.factor(sapply(taxlev,function(x) paste(x[j:j], collapse=";"))))
			V(g1_noedges)$color <- allleg$color
			allleg$Taxa <- NULL
			allleg$Color <-NULL
			leg<-unique(allleg)
			leg <- leg[order(leg$phyla),]
			#order it by name
			l=layout_in_circle(g1,order =order(V(g1)$name))
			#the all group should have the maximum line weight---otherwise something went wrong!
			maxweight = max(E(g1)$weight)
		}
		patientmaxweight = max(E(g1)$weight)
		#add the noedges plus the graph of interest together
		gtp <- g1_noedges + g1
		g1_taxa_df <- taxalist_subbing(paste("/workdir/users/agk85/CDC2/das/",patient,"_taxa.txt",sep=""), j) #filename
                V(gtp)$present <- V(gtp)$taxname %in% g1_taxa_df$taxa		
		#plotting
		#uncolor taxa that are not in it
                V(gtp)$color[!V(gtp)$present]<-NA

		pdf(paste("Org_org_", genetype, "_",patient,"_",contacts, "__",level,"circle.pdf",sep=""),height=20,width=20)
		vsize=4
		#line_scaling=scaling
		plot(gtp,layout=l,vertex.label.cex=2, vertex.label.dist=1, vertex.label=V(gtp)$taxname,vertex.color = V(gtp)$color, vertex.size=vsize,edge.width=E(gtp)$weight_2/line_scaling,edge.color=kewl[i])
		legend("topright",cex=1.75,lty=1,lwd=5,legend=patients,col=c(kewl))
		legend("bottomright",cex = 1.75,pch=21,legend=leg$phyla,pt.bg=c(leg$color))
		legend("bottomleft", cex =3, legend=c("min: 1",paste("half: ", patientmaxweight/2, sep=""),paste("max: ",patientmaxweight, sep="")), lty=1, lwd=c(1,patientmaxweight/line_scaling/2,patientmaxweight/line_scaling) ,bty="n")
		par(new=TRUE)
		dev.off()
	}
	print(level)
	print(densities)	
}




