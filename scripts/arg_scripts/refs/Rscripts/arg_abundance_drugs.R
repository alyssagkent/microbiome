library(ggplot2)
library(RColorBrewer)
library(reshape2)

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/taxonomic_distributions")


#color generation
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 7
kewl<-col_vector[1:(n)]
#pie(rep(1,n), col=kewl)

#import arg rpkm
infile2=read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95.txt",header =F,row.names=1)
genes = as.character(unlist(infile2[1,]))
infile3 = infile2[2:nrow(infile2),]
colnames(infile3)<-genes
df3 = t(infile3)
df4<- as.matrix(as.matrix((df3[,order(colnames(df3))])))
df4_num<-matrix(as.numeric(unlist(df4)),nrow=nrow(df4))
colnames(df4_num)<-colnames(df4)
rownames(df4_num)<-rownames(df4)
df4_withdata <- df4_num[rowSums(df4_num)>0,]
linkerfile = read.csv("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_v_samp_95_95_names_mech.txt",header=T,sep="\t")
df5<-merge(linkerfile,df4_withdata,  by.x="Protein",by.y="row.names", all.y=T)
df5$Protein <-NULL
df5$CARD <-NULL
df5$Resfams <-NULL
df5.melt <- melt(df5, id = c("Cluster", "Name", "Mechanism","Sub_mechanism"))
colnames(df5.melt)<-c("Cluster","Name","Mechanism","Sub_mechanism","Sample","RPKM")
df5.melt$RPKM<-as.numeric(as.character(df5.melt$RPKM))
df5.melt$Patient <- as.vector(sapply(strsplit(as.character(df5.melt$Sample),"-"), function(x) x[[1]]))
df5.melt$TP <- as.vector(sapply(strsplit(as.character(df5.melt$Sample),"-"), function(x) x[[2]]))

test = subset(df5.melt,Sub_mechanism=="Beta-Lactamase")
pdf("ARG_Abundances.pdf", height=20, width=10)
ggplot(test, aes(TP, RPKM, group=Cluster))+
geom_line(aes(color=Name))+
facet_grid(Patient~.,scales="free")
dev.off()

test = subset(df5.melt,Sub_mechanism=="Beta-Lactamase")
pdf("ARG_Abundances_Beta-Lactamase.pdf", height=20, width=10)
ggplot(test, aes(TP, RPKM, group=Cluster))+
geom_line(aes(color=Name))+
facet_grid(Patient~.,scales="free")
dev.off()

van_genes <- c("vanR","vanS","vanT","vanW","vanXYG","vanX","vanY")
test = subset(df5.melt,Name %in% van_genes)
pdf("ARG_Abundances_Vancomycin.pdf", height=20, width=10)
ggplot(test, aes(TP, RPKM, group=Cluster))+
geom_line(aes(color=Name))+
facet_grid(Patient~.,scales="free")
dev.off()
