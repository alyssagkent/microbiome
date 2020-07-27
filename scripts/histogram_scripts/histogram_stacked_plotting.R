#this program will plot the histograms for one sample
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(colorspace)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC2/bins/histogram_figures")
args <- commandArgs(trailingOnly = TRUE)
thresh=args[1]
genetype = args[2] #arg or mge
taxcount = as.numeric(args[3])
inhandle = paste('Together_das_',thresh, '_',genetype,'taxa_together.txt',sep='')
infile_1 = read.csv(inhandle,header=T, sep="\t")
infile_1$Binner=rep("DAS",nrow(infile_1))
#I don't want ot see resident for network
infile1 = subset(infile_1, Thresh=="contacts")

infile = infile1
sample = 'Together'
min_contact = infile_1$Min_contacts[1]
#remap the bintypes
infile$Bintype = factor(infile$Bintype, levels=c("quality","anybin","bin+contig"), labels=c("Quality-bin","Any-bin","Bin or Contig")) 
subfile = subset(infile, Bintype =="Any-bin")
#subfile = subset(infile, Bintype!="Quality-bin")
subfile$Residency = factor(subfile$Residency, levels=c("resident","hic"), labels=c("Resident","Resident or Hi-C"))
s = subset(subfile, TaxaCount>=taxcount & Level == "s__")
s2 = subset(subfile, Level == 's__')

delimname = 's__'
lev = 'species'
df = s
df$TaxaCount <- factor(df$TaxaCount)
df$TaxaCount <- factor(df$TaxaCount, levels=rev(levels(df$TaxaCount)))
plot2<-ggplot(data=df, aes(x=Patient)) +
geom_bar(stat="count", aes(fill=TaxaCount),colour='black',size=0.1)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
scale_fill_manual(values=sequential_hcl(17,"PuBuGn")[min(as.numeric(levels(df$TaxaCount))):max(as.numeric(levels(df$TaxaCount)))])+
facet_grid(.~Residency)+
stat_count(geom="text", aes(label=ifelse(stat(count)>0,stat(count),NA), vjust=-1))


print(paste(genetype,"_ORG_histogram_stacked_species_min1_min", min_contact, "_",taxcount,"+_",sample,".pdf",sep=""))
pdf(paste(genetype,"_ORG_histogram_stacked_species_min1_min", min_contact, "_",taxcount,"+_",sample,".pdf",sep=""), height = 15, width =10)
grid.arrange(plot2, ncol=1, nrow=1, widths=c(1), heights=c(1))
dev.off()

