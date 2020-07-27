library(ggplot2)
library(reshape2)

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.grid = element_blank())}
theme_set(theme_nogrid())
args <- commandArgs(trailingOnly = TRUE)
minreads = args[1]

setwd("/workdir/users/agk85/CDC2/args/patric_comparisons")
topargs= read.csv("/workdir/users/agk85/CDC2/args/patric_comparisons/all_clusters_topargs_mergednames.txt", sep="\t",header = T)
topargs$Mergedname <- NULL

arginfo= read.csv("/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech.txt", sep="\t",header = T)
patric <- read.csv("/workdir/users/agk85/CDC2/args/patric_comparisons/Patric_arg_taxonomies.txt",header=F, sep="\t")
hic <- read.csv(paste("/workdir/users/agk85/CDC2/args/patric_comparisons/Hic_arg_", minreads, "_taxonomies.txt",sep=""),header=F, sep="\t")
base <- read.csv(paste("/workdir/users/agk85/CDC2/args/patric_comparisons/Base_arg_", minreads, "_taxonomies.txt",sep=""),header=F, sep="\t")

#combine base and hic # might eventually change this
hicbase <- unique(rbind(hic,base))
hicbase$Origin <- rep("Study",nrow(hicbase))

patric$Origin <- rep("PATRIC",nrow(patric))
df <- rbind(hicbase,patric)
header = c("Cluster","Taxonomy","Type")
colnames(df)<- header

#add on the arg info
dfai <- merge(df,arginfo,by=c("Cluster"))

#keep only the topargs
topdf <- merge(topargs,dfai,by=c("Cluster"))

levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
delims = c("k__","p__","c__","o__","f__","g__","s__")
### make all of the Taxonomy all split up
for (i in seq(7)){
#get the level name
level = levels[i]
#get the delimiter
delim = delims[i]
#get the stuff before the delimiter and after and the actual level taxa, then merge front and taxa
topdf$Taxafront <- sapply(strsplit(as.character(topdf$Taxonomy),delim), function(x) x[[1]])
topdf$Taxaback <- sapply(strsplit(as.character(topdf$Taxonomy),delim), function(x) x[[2]])
topdf$Taxa <- sapply(strsplit(as.character(topdf$Taxaback),";"), function(x) x[[1]])
topdf$Taxafull <- paste(topdf$Taxafront, topdf$Taxa)
top <- topdf
top$Taxonomy <- NULL
# top$Taxafull <- NULL
topdfuniq<-unique(top)
topdffilled <- subset(topdfuniq, Taxa != "")
topdffilled <- topdfuniq
outhandle = paste("ARG_Topargs_",minreads,"_",level, "_Distributions.pdf",sep="")
pdf(outhandle,height=2.5*i,width=i*3.5)
g<-ggplot(topdffilled, aes(x=Taxafull))+
geom_bar(stat="count",aes(fill=Alltypes))+
facet_grid(Type~.,scales="free_y")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
scale_fill_manual(values=c("brown1","burlywood1","cadetblue1","navy","plum4"))
print(g)
dev.off()
}


