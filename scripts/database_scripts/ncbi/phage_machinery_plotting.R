#plotting phage-ncbi distribution
library(ggplot2)
library(reshape2)

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.grid = element_blank())}
theme_set(theme_nogrid())

setwd("/workdir/users/agk85/CDC2/bins/patric_figures")

args <- commandArgs(trailingOnly = TRUE)
minreads = args[1]
genetype = args[2]
support = args[3]
reference = args[4]

infile = read.table(paste(genetype, "_", minreads, "_",support, "_",reference, "_distribution.txt",sep=""),sep="\t", header=T)
infile$nohit<-NULL
infile$total<-NULL
infile$overlap<-NULL
infile$no_overlap <- NULL
infile$expand <- NULL
df.melt <- melt(infile, id=c("delimname"))

colnames(df.melt)<-c("Level","Category","Count")
df.melt$Level <- factor(df.melt$Level, levels=c('species','genus','family','order','class','phylum','kingdom'))
df.melt$Category <- factor(df.melt$Category, levels=c("unannotated_level",'num1','num2',"num4","num5_2","num7","num9","num3","num5","num6", "num8","num8_2"))

pdf(paste("phage-machinery_", minreads, "_", support, "ncbi_distribution.pdf",sep=""), height = 5, width=8)
ggplot(df.melt, aes(x=Level, y=Count))+
geom_bar(aes(fill=Category),stat="identity")+
ylab("Phage Machinery Gene Cluster Count")+
xlab("Taxonomic Level")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
scale_fill_manual(values=c("unannotated_level"="gray65","num1"="gray80","num2"="gray95",
"num3"="deepskyblue","num5"="blue","num6"="blue2","num8" = "green","num8_2"="navy",
"num4"="darkorange","num5_2"="darkorange2","num7"="darkorange3", "num9"="darkorange4"), 
labels = c("num1"="No NCBI Hits or Viral-One HiC Taxon",
"num2"="No NCBI Hits or Viral-Multiple HiC Taxa",
"unannotated_level"="Unannotated NCBI or HiC Taxa at this Level",
"num3"="One NCBI Taxon-One HiC Taxon-Overlapping",
"num5"="One NCBI Taxon-Multiple HiC Taxa-Overlapping",
"num6"="Multiple NCBI Taxa-One HiC Taxon-Overlapping",
"num4"="One NCBI Taxon-One HiC Taxon-Not Overlapping",
"num5_2"="One NCBI Taxon-Multiple HiC Taxa-Not Overlapping",
"num7"="Multiple NCBI Taxa-One HiC Taxon-Not Overlapping",
"num8"="Multiple NCBI Taxa-Multiple HiC Taxa-Overlapping Completely",
"num8_2"="Multiple NCBI Taxa-Multiple HiC Taxa-Overlapping",
"num9"="Multiple NCBI Taxa-Multiple HiC Taxa-Not Overlapping"))+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#1 There was no known hit., we assign 1 taxa. (putative host-specific phage)
#2 There was no known hit, we assign multiple taxa (putative broad-range phage)

#3 There was 1 hit, we get the same (confirms host-specific phage)
#4 There was 1 hit, we get 1 different hit (putative broad range, unknown before, or crappy hit)
#5 There was 1 hit, we get more hits, including the BLAST hit (putative broad range, unknown before, or crappy hit)
#5_2 There was 1 hit, we get more hits, not including the BLAST hit (putative broad range, unknown before, or crappy hit)

#6 There was multiple hits, we get some overlap (confirms at least a bit of what we know)
#7 There was multiple hits, we get 1 or more, none of which overlap (expand our knowledge of organisms it can invade).

