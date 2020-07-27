#Hi-C trans interactions mapped in igraph
#
library(igraph)
library(qgraph)
library(ggplot2)
setwd('/workdir/users/agk85/CDC/newhic/plotting/')
workdir='/workdir/users/agk85/CDC'
args <- commandArgs(trailingOnly = TRUE)
identity = args[1]

df <- data.frame(Name=character(),
	Type=character(),
	Group=character(),
	Percents=numeric(),
	MGE_counts=numeric(),
	MGE_totals=numeric(),
	Cutoff=numeric()) 

argorgdf <- data.frame(Name=character(),
	v1 = character(),
	v2 = character(),
	weight = numeric())

names = read.table(paste(workdir,'/HicDesign.txt',sep=""), header=F)
names = as.character(names$V1)
for (name in names){
#trans interactions from hic data
#noeuks
data <- read.table(paste(workdir, "/newhic/mapping/", name,"/",name, "_trans_primary_", identity, "_noeuks.txt",sep=""), sep="\t")
colnames(data)<-c("ID1","Flag1","Contig1","Start1","Mapq1","CIGAR1","ID2","Flag2","Contig2","Start2","Mapq2","CIGAR2")

annotations <- read.table(paste(workdir, "/combo_tables/metagenomes4/", name,"_master_scf_table_binary.txt",sep=""),sep="\t",header=T,row.names=1)

amph <- read.table(paste(workdir, "/combo_tables/metagenomes4/", name,"_master_scf_table.txt",sep=""),quote = "", sep="\t",header=T,row.names=1)

plasmids <- annotations$Plasmid_Finder + annotations$Full_Plasmids + annotations$Relaxase + annotations$Aclame_plasmid + annotations$Imme_plasmid + annotations$Plasmid_pfams + annotations$Pfams_plasmid
args <- annotations$Resfams + annotations$Perfect + annotations$CARD
is <- annotations$ISEScan + annotations$Imme_transposon + annotations$Pfams_transposon
phage <- annotations$Phage_Finder  + annotations$Phaster + annotations$Imme_phage + annotations$Aclame_phage + annotations$Virus + annotations$Pfams_phage
org <- annotations$Amphora + annotations$Rnammer + annotations$Campbell + annotations$Metaphlan + annotations$GSMer
#TODO make sure this little bit works....do I have to do something??
orgname<-amph$Amphora.best_taxonomy.confidence..
argname<-paste(amph$Resfams, amph$Perfect, amph$CARD, sep="")
annot<- data.frame(plasmids, args, is, phage, org,orgname,argname )
rownames(annot)<-rownames(annotations)

#testing purposes a subset of the interactions
data_s <- data
#get the links c1 is first contig and c2 is second contig
a <- data.frame(as.character(data_s$Contig1),as.character(data_s$Contig2))
g<-graph_from_data_frame(a, directed=FALSE)
#give it weights for the edges
E(g)$weight<-1

#simplify by collapsing multiple edges and summing them
g_simp <- simplify( g, remove.multiple=T, edge.attr.comb=c(weight="sum") )

#normalize the width to the biggest width (so you don't have whoppers)
E(g_simp)$width <- E(g_simp)$weight   #E(g_simp)$weight/(max(E(g_simp)$weight)/30)

#set a layout 
#l <- layout_with_kk(g_simp)
#plot(g_simp,vertex.label=NA,vertex.size=5,vertex.color=c( "white", "red")[1+(V(g_simp)$plasmid==1)] )

#l <- layout_with_kk(g_simp)
#plot(g_simp, layout=l, vertex.label=NA, vertex.size=5)
hist(E(g_simp)$weight,breaks=35)
mean(E(g_simp)$weight)


##############
cutoffs=c(2)
for (cut.off in cutoffs){
g_one <- delete_edges(g_simp, E(g_simp)[weight<cut.off])
V(g_one)$deg <- degree(g_one, mode="all")                                              
g_one=delete.vertices(g_one,which(degree(g_one)<1))                              

#make sure sort is False because when  you just tack it on the vertex attributes it has tobe the same order
b<- merge(as.data.frame(V(g_one)$name), annot,by.x="V(g_one)$name",by.y=0,all.x=T,sort=F)
dim(b)
dim(data)

rownames(b)<-b$`V(g_one)$name`
b$`V(g_one)$name` <- NULL
b_orgname <- b$orgname
b$orgname <- NULL

b_argname <- b$argname
b$argname <- NULL

b$sums <- rowSums(b)
V(g_one)$plasmid<-b$plasmids
V(g_one)$phage<-b$phage
V(g_one)$arg<-b$args
V(g_one)$is<-b$is
V(g_one)$org<-b$org
V(g_one)$deg <- degree(g_one, mode="all")
V(g_one)$sums <- rowSums(b)
V(g_one)$orgname <- as.character(b_orgname)
V(g_one)$argname <- as.character(b_argname)
g_one_annot=delete_vertices(g_one,which(V(g_one)$sums<1))
V(g_one_annot)$deg <- degree(g_one_annot, mode="all")
g_one_annot=delete.vertices(g_one_annot,which(degree(g_one_annot)<1))

#########################################################################################
#get the percentage captured 
mge_percents = rep(0,12)
mge_captured = rep(0,12)
mge_total = rep(0,12)
#percent mgm mge's captured in trans interactions 
mge_percents[1]<-100*sum(V(g_one)$plasmid>0)/sum(plasmids>0)
mge_percents[2]<-100*sum(V(g_one)$phage>0)/sum(phage>0)
mge_percents[3]<-100*sum(V(g_one)$is>0)/sum(is>0)
mge_percents[4]<-100*sum(V(g_one)$arg>0)/sum(args>0)
mge_percents[5]<-100*sum(V(g_one)$org>0)/sum(org>0)
mge_percents[6]<-100*vcount(g_one)/length(annotations$Plasmid_Finder)
mge_captured[1]<-sum(V(g_one)$plasmid>0)
mge_captured[2]<-sum(V(g_one)$phage>0)
mge_captured[3]<-sum(V(g_one)$is>0)
mge_captured[4]<-sum(V(g_one)$arg>0)
mge_captured[5]<-sum(V(g_one)$org>0)
mge_captured[6]<-vcount(g_one)
mge_total[1]<-sum(plasmids>0)
mge_total[2]<-sum(phage>0)
mge_total[3]<-sum(is>0)
mge_total[4]<-sum(args>0)
mge_total[5]<-sum(org>0)
mge_total[6]<-length(annotations$Plasmid_Finder)
#percent mgm mge's captured in trans interactions connected to 
mge_percents[7]<-100*sum(V(g_one_annot)$plasmid>0)/sum(plasmids>0)
mge_percents[8]<-100*sum(V(g_one_annot)$phage>0)/sum(phage>0)
mge_percents[9]<-100*sum(V(g_one_annot)$is>0)/sum(is>0)
mge_percents[10]<-100*sum(V(g_one_annot)$arg>0)/sum(args>0)
mge_percents[11]<-100*sum(V(g_one_annot)$org>0)/sum(org>0)
mge_percents[12]<-100*vcount(g_one_annot)/length(annotations$Plasmid_Finder)
mge_captured[7]<-sum(V(g_one_annot)$plasmid>0)
mge_captured[8]<-sum(V(g_one_annot)$phage>0)
mge_captured[9]<-sum(V(g_one_annot)$is>0)
mge_captured[10]<-sum(V(g_one_annot)$arg>0)
mge_captured[11]<-sum(V(g_one_annot)$org>0)
mge_captured[12]<-vcount(g_one_annot)
mge_total[7]<-sum(plasmids>0)
mge_total[8]<-sum(phage>0)
mge_total[9]<-sum(is>0)
mge_total[10]<-sum(args>0)
mge_total[11]<-sum(org>0)
mge_total[12]<-length(annotations$Plasmid_Finder)

sampleids <- rep(name,12)
cfs <- rep(cut.off, 12)
type <- c("plasmid","phage","is","arg","org","trans","plasmid","phage","is","arg","org","trans")
grp <-c(rep("Captured.by.trans.interactions",6),rep("Both.annotated",6))
newdf<-data.frame(sampleids, type, grp,mge_percents,mge_captured,mge_total, cfs)
df <- rbind(df,newdf)

#####################################################
#Plotting Bubble Plot
#####################################################
clean_argname<-sapply(strsplit(V(g_one_annot)$argname,"\\|"), `[`, 2)
cleaner_argname<-sapply(strsplit(clean_argname, "ARG_"),`[`,1)
cleaner_argname[is.na(cleaner_argname)]<-""


clean_orgname<-sapply(strsplit(V(g_one_annot)$orgname,"\\_"), `[`, 2)
clean_orgname[is.na(clean_orgname)]<-""

#l <- layout_nicely(g_one_annot)
#l <- layout_with_fr(g_one_annot)
#plot(g_one_annot,layout=l,vertex.label=NA,vertex.size=7,vertex.color=c( rgb(0, 0, 0, alpha = 0), "red")[1+(V(g_one_annot)$plasmid>0)] )
g_edge <- g_one_annot
E(g_edge)$weight<-1
graph.to.plot <- g_edge
l <- layout_with_mds(graph.to.plot)


########change this back when you want to replot data 
pdf(paste("HIC_trans_nicely_annot_withlink_", name,"_", identity, "_",cut.off, ".pdf",sep=""),height=200, width=200)
 
 plot(graph.to.plot,layout=l,edge.color=NA,vertex.label=NA,vertex.size=5.5,vertex.color=c(NA, "red")[1+(V(graph.to.plot)$plasmid>0)] )
 par(new=TRUE)
 plot(graph.to.plot,layout=l,vertex.frame.color=NA,edge.color=NA,vertex.label=NA,vertex.size=5,vertex.color=c(NA, "blue")[1+(V(graph.to.plot)$phage>0)] )
 par(new=TRUE)
 plot(graph.to.plot,layout=l,vertex.frame.color=NA,edge.color=NA,vertex.label=NA,vertex.size=4.5,vertex.color=c(NA, "green")[1+(V(graph.to.plot)$is>0)] )
 par(new=TRUE)
 plot(graph.to.plot,layout=l,vertex.frame.color=NA,edge.color=NA,vertex.label=NA,vertex.size=4,vertex.color=c(NA,"purple")[1+(V(graph.to.plot)$org>0)] )
 par(new=TRUE)
 plot(graph.to.plot,layout=l,vertex.frame.color=NA,vertex.size=3,vertex.label=paste(cleaner_argname,clean_orgname,sep=""),mark.border=NA,
vertex.color=c(NA, "orange")[1+(V(graph.to.plot)$arg>0)],edge.width = E(graph.to.plot)$weight )
legend("bottomright",c("Plasmid","Phage","Transposon","Lineage-specific gene","AR gene"), col=c("red","blue","green","purple","orange"))
dev.off()

#################arg to organism ##################################

a <-E(g_one_annot)[inc(V(g_one_annot)[org>0])]
g2 <- subgraph.edges(g_one_annot, a)
b <-E(g2)[inc(V(g2)[arg>0])]
g3<- subgraph.edges(g2, b)

###how to get interactions out!??
V(g3)$name<-paste(V(g3)$name, V(g3)$argname, V(g3)$orgname,sep="#")
argorgdat <-data.frame(as_edgelist(g3))
argorgdat$weight <- E(g3)$weight
samp <- rep(name, dim(argorgdat)[1])
argorgdat$samp <- samp
#argorgdf <- rbind(argorgdf,argorgdat)

###########using the argorg data
#so remove things that don't have an annotation
temp <- argorgdat
temp$X1_value <- grepl("amphora", temp$X1) + grepl("ARG", temp$X1)
temp$X2_value <- grepl("amphora", temp$X2) + grepl("ARG", temp$X2)
temp$X1[temp$X1_value<1]<-NA
temp$X2[temp$X2_value<1]<-NA
t<-temp[complete.cases(temp), ]
tnew <-data.frame(t$samp,t$weight, t$X1, t$X2)
colnames(tnew)<-c("Sample","Links","Contig1","Contig2")
write.csv(tnew, paste("Org_ARG_links_", name, "_", identity, "_", cut.off, ".csv", sep=""))

}
}

write.csv(df, paste("CDC_percent_capture_", identity, "id.csv",sep=""))
