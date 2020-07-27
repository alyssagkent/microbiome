#this program will plot the histograms for one sample
library(ggplot2)
library(gridExtra)
library(reshape2)
library(reshape)

theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC2/bins")
args <- commandArgs(trailingOnly = TRUE)
thresh=args[1]
genetype = args[2] #arg or mge

inhandle = paste('all_bin_cleanliness_all.txt',sep='')
df = read.csv(inhandle,header=T, sep="\t")

#proportions
df$intra_bin_cis

df.melt<-melt(df, by=


plot2<-ggplot(data=df, aes(TaxaCount)) +
geom_histogram(binwidth = 1, size=1)+
labs(title=paste(sample))+
xlab(paste("Number of ", lev)) +
ylab(paste(genetype, " contig Taxa Counts",sep=""))+
scale_x_continuous(breaks=seq(2,max(df$TaxaCount),1))+
facet_grid(Patient~Residency)+
#facet_grid(Patient~Binner+Thresh+Bintype+Residency)+
stat_bin(geom="text", aes(label=ifelse(stat(count)>0,stat(count),NA), vjust=-1))

pdf(paste("Cleanliness_bin_numbers_", mge, "_",".pdf",sep=""), height = 15, width =10)
grid.arrange(plot2, ncol=1, nrow=1, widths=c(1), heights=c(1))
dev.off()

