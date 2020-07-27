#this program will plot the histograms for one sample
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(colorspace)
library(reshape2)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank()
    )
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC2/bins/timelapse")
args <- commandArgs(trailingOnly = TRUE)
thresh=args[1]
genetype = args[2] #arg or mge
taxlevel = args[3]

inhandle = paste('timelapse_', genetype, '_org_',taxlevel,'_',thresh,'.txt',sep='')
infile = read.csv(inhandle,header=T, sep="\t")
infile$Total_T1 =NULL
infile$Total_T2 = NULL
infile$Uniqelement1 = NULL
infile$Uniqelement2 = NULL
infile$Level = NULL

subfile <- subset(infile, Timepoint_of_comparison==1)
subfile$Timepoint_of_comparison <- NULL
df = melt(subfile, id.vars=c("T1","T2"))

pdf(paste("Timelapse_",genetype, '_org_', taxlevel, '_',thresh,'_lines.pdf',sep=''),height=7,width=10)
ggplot(df, aes(x=T2, y=value, group =variable))+
geom_line(aes(color=variable))+
geom_point()+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
scale_color_manual(values=c("#800080","#FF0000","#0000CD"))+
facet_wrap(.~T1,scales="free_x",nrow=3)
dev.off()

