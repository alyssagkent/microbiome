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
setwd("/workdir/users/agk85/CDC2/resampling")
args <- commandArgs(trailingOnly = TRUE)
inhandle=args[1]
outhandle=args[2]
infile = read.table(inhandle, header=T, sep="\t")

pdf(outhandle,height=6,width=8)
ggplot(infile,aes(x=reads, y=unique_connections,group=sample,colour=patient))+
	geom_line()+
	geom_point()+
	xlab("Reads")+
	ylab("Unique connections")
dev.off()


