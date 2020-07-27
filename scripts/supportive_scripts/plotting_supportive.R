library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())

setwd("/workdir/users/agk85/CDC2/bins_hicsupport")


#for each patient show the number of connections for each category facet_grid

header=c("minread","patient","thresh","bintype","residency","level","geneid","gene_mech","taxa","genetype")
infile1 = read.table("Together_das_2_argtaxa_linelist.txt",header=F,sep="\t")
infile2 = read.table("Together_das_5_argtaxa_linelist.txt",header=F,sep="\t")
infile3 = read.table("Together_das_2_mgetaxa_linelist.txt",header=F,sep="\t")
infile4 = read.table("Together_das_5_mgetaxa_linelist.txt",header=F,sep="\t")

g1=rep("arg",nrow(infile1))
g2=rep("arg",nrow(infile2))
g3=rep("mge",nrow(infile3))
g4=rep("mge",nrow(infile4))

infile1$genetype=g1
infile2$genetype=g2
infile3$genetype=g3
infile4$genetype=g4

df<- rbind(infile1,infile2,infile3,infile4)
colnames(df)<-header
levs = c("s__","g__")
for (lev in levs){
dfsub <- subset(df, thresh=="contacts" & bintype=="anybin" & level==lev)

pdf(paste("hic_support_",lev, ".pdf",sep=""),height=10, width=10)
g <- ggplot(dfsub, aes(x=patient))+
geom_bar(stat="count",aes(fill=residency),colour='black',size=0.1)+
theme(axis.text.x=element_text(angle=90,hjust=1))+
scale_fill_brewer(palette = "Set3")+
facet_grid(genetype~minread,scales="free_y")+
stat_count(geom="text",aes(label=ifelse(stat(count)>0,stat(count),NA),vjust=-0.3))
print(g)
dev.off()
a = dfsub %>% group_by(patient,minread,genetype,residency) %>% tally()
write.table(a,paste("hic_support_numbers_",lev, ".txt",sep=""),sep="\t",row.names=F,quote=F)
}

