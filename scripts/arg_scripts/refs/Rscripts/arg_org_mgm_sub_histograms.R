#plot histograms for arg org and org arg relationships
#updated to only include things that have taxonomy down to species to be counted
library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(scales)
library(circlize)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mgm_sub_histograms")

patients = c("B314_","B316_","B320_","B331_","B335_","B357_","B370_","")
#patients = c("")
delims = c('species','genus','family','order','class','phylum','kingdom')
figs <- vector("list", 24) 
keepcount <- vector("list", 7)
stats <- vector("list", 8)
for (i in seq(length(delims))){
figs <- vector("list", 24)
for (j in seq(length(patients))){
patient = patients[j]
print(patient)
delimname = delims[i]
print(delimname)
indata0 = read.csv(paste("arg_orgcounts_", patient,delimname, "0.txt",sep=""), sep="\t",  header=F)
indata1 = read.csv(paste("arg_orgcounts_",patient,delimname, "1.txt",sep=""), sep="\t", row.names=1, header=F)
indata2 = read.csv(paste("arg_orgcounts_",patient,delimname, "2.txt",sep=""), sep="\t", row.names=1, header=F)
combocounts = cbind(indata0, indata1, indata2)
colnames(combocounts)<-c("Cluster","arg0","count0","arg1","count1","arg2","count2")
#assess mean sd median of counts
mean_all= c(mean(combocounts$count0), mean(combocounts$count1),mean(combocounts$count2))
median_all=c(median(combocounts$count0), median(combocounts$count1),median(combocounts$count2))
sd_all = c(sd(combocounts$count0), sd(combocounts$count1),sd(combocounts$count2))

combocounts[combocounts == 0] <- NA
colnames(combocounts)<-c("Cluster","arg0","count0","arg1","count1","arg2","count2")
combo<-combocounts
keepcount[[j]] <- combo
notzero0 <- round(100*sum(!(is.na(combo$count0)))/nrow(combo),2)
notzero1 <- round(100*sum(!(is.na(combo$count1)))/nrow(combo),2)
notzero2 <- round(100*sum(!(is.na(combo$count2)))/nrow(combo),2)


mean_nozero= c(mean(combocounts$count0,na.rm=T), mean(combocounts$count1,na.rm=T),mean(combocounts$count2,na.rm=T))
median_nozero=c(median(combocounts$count0,na.rm=T), median(combocounts$count1,na.rm=T),median(combocounts$count2,na.rm=T))
sd_nozero = c(sd(combocounts$count0,na.rm=T), sd(combocounts$count1,na.rm=T),sd(combocounts$count2,na.rm=T))
max_nozero = c(max(combocounts$count0,na.rm=T), max(combocounts$count1,na.rm=T),max(combocounts$count2,na.rm=T))
print(max_nozero)
statvec = c(mean_all, median_all, sd_all, mean_nozero, median_nozero, sd_nozero)
stats[[j]] = statvec
#determine shared ylimit
temp0<-ggplot(data=combo, aes(count0)) + 
geom_histogram(binwidth = 1)
temp1<-ggplot(data=combo, aes(count1)) + 
geom_histogram(binwidth = 1)
temp2<-ggplot(data=combo, aes(count2)) + 
geom_histogram(binwidth = 1)
k0 = ggplot_build(temp0)
k1 = ggplot_build(temp1)
k2 = ggplot_build(temp2)

m0 = max(k0$data[[1]][2])
m1= max(k1$data[[1]][2])
m2= max(k2$data[[1]][2])

m = max(m0, m1, m2)

a<-ggplot(data=combo, aes(count0)) + 
geom_histogram(binwidth = 1, size=2)+
ylim(0,m+25)+
xlim(0.5,5)+
labs(title=paste(patient, " Assembly only"))+
xlab(paste("Number of ", delimname)) +
ylab("ARG frequency")+
annotate("text", label = paste(notzero0, "%"), x = max(combo$count0,na.rm=T), y = .25*m , size = 3)

b<-ggplot(data=combo, aes(count1)) + 
geom_histogram(binwidth = 1)+
ylim(0,m+25)+
xlim(0.5,26)+
labs(title=paste("Assembly + 1 link"))+
xlab(paste("Number of ", delimname)) +
ylab("ARG frequency")+
annotate("text", label = paste(notzero1, "%"), x = max(combo$count1,na.rm=T), y = .25*m , size = 3)


c<-ggplot(data=combo, aes(count2)) + 
geom_histogram(binwidth = 1)+
ylim(0,m+25)+
xlim(0.5,61)+
labs(title=paste(patient, " Assembly + 2 links"))+
xlab(paste("Number of ", delimname)) +
ylab("ARG frequency")+
annotate("text", label = paste(notzero2, "%"), x = max(combo$count2,na.rm=T), y = .25*m , size = 3)

figs[[j]]= a
figs[[j+8]]= b
figs[[j+16]]= c

}

pdf(paste("All_ARG_ORG_MGM_SUB_histogram_no0_",delimname,"_with2.pdf",sep=""), height = 20, width =15)
grid.arrange(figs[[1]], figs[[9]], figs[[17]],
	figs[[2]], figs[[10]],figs[[18]],
	figs[[3]], figs[[11]],figs[[19]],
	figs[[4]], figs[[12]],figs[[20]],
	figs[[5]], figs[[13]],figs[[21]],
	figs[[6]], figs[[14]],figs[[22]],
	figs[[7]], figs[[15]],figs[[23]],
	figs[[8]], figs[[16]],figs[[24]],
	ncol=3, nrow=8, widths=c(2.5, 5, 8), heights=c(rep(2,8)))
dev.off()

}

