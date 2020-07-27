#this program will plot the histograms for one sample
library(ggplot2)
library(gridExtra)

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
plot2<-ggplot(data=df, aes(TaxaCount)) +
geom_histogram(binwidth = 1)+
labs(title=paste(sample))+
xlab(paste("Number of ", lev)) +
ylab(paste(genetype, " contig Taxa Counts",sep=""))+
scale_x_continuous(breaks=seq(taxcount,max(df$TaxaCount)))+
facet_grid(Patient~Residency)+
#facet_grid(Patient~Binner+Thresh+Bintype+Residency)+
stat_bin(geom="text", aes(label=ifelse(stat(count)>0,stat(count),NA), vjust=-1))

print(paste(genetype,"_ORG_histogram_species_min1_min", min_contact, "_",taxcount,"+_",sample,".pdf",sep=""))
pdf(paste(genetype,"_ORG_histogram_species_min1_min", min_contact, "_",taxcount,"+_",sample,".pdf",sep=""), height = 15, width =10)
grid.arrange(plot2, ncol=1, nrow=1, widths=c(1), heights=c(1))
dev.off()

df$Counter <- rep(1,nrow(df))
df$Taxaset <- NULL
agg <-aggregate(Counter ~ Patient+Thresh+Bintype+Residency+Level+Binner, df, sum)
agg$Min_contact <- rep(min_contact, nrow(agg))
agg$Total_args <- rep(df$Total_args[1],nrow(agg))
print(agg)
write.table(agg, paste(genetype,"_Numbers_","min1_min",min_contact, "_",sample,"_2+.txt",sep=""), row.names=F,quote=FALSE,sep="\t")


agg <-aggregate(Counter ~ Patient+Thresh+Bintype+Residency+Level+Binner+factor(TaxaCount), df, sum)
agg$Min_contact <- rep(min_contact, nrow(agg))
agg$Total_args <- rep(df$Total_args[1],nrow(agg))
print(agg)
write.table(agg, paste(genetype, "_Taxa_Count_Numbers_","min1_min",min_contact, "_",sample,"_2+.txt",sep=""),row.names=F, quote=FALSE,sep="\t")

agg_sum <-aggregate(TaxaCount ~ Patient+Thresh+Bintype+Residency+Level+Binner, df, sum)
agg_sum$metric <- rep("sum",nrow(agg_sum))
agg_mean <-aggregate(TaxaCount ~ Patient+Thresh+Bintype+Residency+Level+Binner, df, mean)
agg_mean$metric <- rep("mean",nrow(agg_mean))
agg_median<-aggregate(TaxaCount ~ Patient+Thresh+Bintype+Residency+Level+Binner, df, median)
agg_median$metric <- rep("median",nrow(agg_median))
agg <- rbind(agg_sum,agg_mean, agg_median)
write.table(agg, paste(genetype,"_Metrics_","min1_min",min_contact, "_",sample,"_2+.txt",sep=""), row.names=F,quote=FALSE,sep="\t")


s2$Counter <- rep(1,nrow(s2))
s2$Taxaset <- NULL
agg <-aggregate(Counter ~ Patient+Thresh+Bintype+Residency+Level+Binner, s2, sum)
agg$Min_contact <- rep(min_contact, nrow(agg))
agg$Total_args <- rep(s2$Total_args[1],nrow(agg))
print(agg)
write.table(agg, paste(genetype,"_Numbers_","min1_min",min_contact, "_",sample,"_1+.txt",sep=""), row.names=F,quote=FALSE,sep="\t")


agg <-aggregate(Counter ~ Patient+Thresh+Bintype+Residency+Level+Binner+factor(TaxaCount), s2, sum)
agg$Min_contact <- rep(min_contact, nrow(agg))
agg$Total_args <- rep(s2$Total_args[1],nrow(agg))
print(agg)
write.table(agg, paste(genetype, "_Taxa_Count_Numbers_","min1_min",min_contact, "_",sample,"_1+.txt",sep=""),row.names=F, quote=FALSE,sep="\t")

agg_sum <-aggregate(TaxaCount ~ Patient+Thresh+Bintype+Residency+Level+Binner, s2, sum)
agg_sum$metric <- rep("sum",nrow(agg_sum))
agg_mean <-aggregate(TaxaCount ~ Patient+Thresh+Bintype+Residency+Level+Binner, s2, mean)
agg_mean$metric <- rep("mean",nrow(agg_mean))
agg_median<-aggregate(TaxaCount ~ Patient+Thresh+Bintype+Residency+Level+Binner, s2, median)
agg_median$metric <- rep("median",nrow(agg_median))
agg <- rbind(agg_sum,agg_mean, agg_median)
write.table(agg, paste(genetype,"_Metrics_","min1_min",min_contact, "_",sample,"_1+.txt",sep=""), row.names=F,quote=FALSE,sep="\t")

