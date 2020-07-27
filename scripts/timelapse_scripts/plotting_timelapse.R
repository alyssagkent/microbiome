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

df = melt(infile, id.vars=c("T1","T2","Level","Timepoint_of_comparison"))
pdf(paste("Timelapse_",genetype, '_org_', taxlevel, '_',thresh,'.pdf',sep=''),height=12,width=30)
ggplot(df, aes(x=T2, y=value,group=variable,fill=variable))+
 geom_bar(stat="identity",position = "dodge")+
 theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 scale_fill_manual(values=c("#800080","#F08080","#6495ED","#FF0000","#0000CD"))+
 facet_grid(Timepoint_of_comparison~T1, scales = "free_x")
dev.off()

#
#
#"#800000", "#191970",

#eventual stacking
#ggplot(df.m, aes(strain), ylim(-500:500)) + 
#geom_bar(data = subset(df.m, variable == "count.up"), 
#   aes(y = value, fill = condition), stat = "identity", position = "dodge") +
#geom_bar(data = subset(df.m, variable == "count.down"), 
#   aes(y = -value, fill = condition), stat = "identity", position = "dodge") + 
#geom_hline(yintercept = 0,colour = "grey90")

