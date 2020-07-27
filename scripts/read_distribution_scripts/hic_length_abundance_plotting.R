library(ggplot2)
library(gplots)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
library(viridis)
library(formattable)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

setwd("/workdir/users/agk85/CDC2/read_distributions/")
indata = read.csv("All_length_hichits_mobility_abundance.txt",header = T,sep="\t")
#colnames(indata)<- c("Scfid","Sample","Length","Hic_hits","Mobile")
indata$Mobility <- as.character(indata$Mobility)
#subdata <- subset(indata, Sample=="B314-1")
colnames(indata)<-c("Scfid","Name", "Length","Hic_hits", "Mobility", "RPKM", "Coverage", "Cis_hic_hits")
subdata <- subset(indata, Coverage>80 )

#goodsamps <- c("B335-1","B335-2","B335-3", "B357-1","B357-2","B357-3","B357-4","B357-5","B357-6","B370-1","B370-2","B370-3")
#goodsamps <- c("B335-1","B335-2","B335-3")
#goodsamps <- c("B370-1","B370-2","B370-3")
#goodsamps <- c("B357-1","B357-2","B357-3","B357-4","B357-5","B357-6")
#subdata <- subset(subdata, Name %in% goodsamps)


r1 <- subset(subdata, Mobility == "1")
r2 <- subset(subdata, Mobility == "0")

cor.test(r1$Length, r1$Hic_hits)
cor.test(r1$RPKM, r1$Hic_hits)
cor.test(r2$Length, r2$Hic_hits)
cor.test(r2$RPKM, r2$Hic_hits)
cor.test(subdata$Length, subdata$Hic_hits)
cor.test(subdata$RPKM, subdata$Hic_hits)

#density plot of each

sub2 <- subset(subdata, Mobility==0 & Hic_hits>0)
dense2 <- ggplot(sub2, aes(RPKM, Hic_hits))+
	geom_bin2d()+
	scale_x_log10(labels=comma)+
	scale_y_log10(labels=comma)+
	ggtitle("Non-Mobile Abundance vs Trans HiC")+
	scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

sub3 <- subset(subdata, Mobility==1& Hic_hits>0)
dense3 <- ggplot(sub3, aes(RPKM, Hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
	ggtitle("Mobile Abundance vs Trans HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

#density plot of each
dense5 <- ggplot(sub2, aes(Length, Hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
	ggtitle("Non-Mobile Length vs Trans HiC")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense6 <- ggplot(sub3, aes(Length, Hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
	ggtitle("Mobile Length vs Trans HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")


sub5<- rbind(sub2,sub3)
plot_top_rpkm <- ggplot(sub5, aes(RPKM, fill=Mobility)) +
  	geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
  	ggtitle("Non-Mobile vs Mobile Abundance vs Trans HiC")+
	scale_fill_manual(values = c("orange", "purple")) +
  	theme(legend.position = "none")

plot_top_length <- ggplot(sub5, aes(Length, fill=Mobility)) +
  geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
	ggtitle("Non-Mobile vs Mobile Length vs Trans HiC")+
  	scale_fill_manual(values = c("orange", "purple")) +
  	theme(legend.position = "none")


dense8 <- ggplot(sub2, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense9 <- ggplot(sub3, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

pdf("Densities_scaffolds_withtrans.pdf",height=12,width=10)
grid.arrange(dense2, dense5,dense3, dense6, plot_top_rpkm, plot_top_length, nrow = 3, ncol=2)
dev.off()

pdf("Length_vs_Abundance_scaffolds_withtrans.pdf",height=12,width=5)
grid.arrange(dense8,dense9, nrow=3,ncol=1)
dev.off()
####################################################################################with all scaffolds

sub2 <- subset(subdata, Mobility==0)
dense2 <- ggplot(sub2, aes(RPKM, Hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Abundance vs Trans HiC")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

sub3 <- subset(subdata, Mobility==1)
dense3 <- ggplot(sub3, aes(RPKM, Hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Abundance vs Trans HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

#density plot of each
dense5 <- ggplot(sub2, aes(Length, Hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Length vs Trans HiC")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense6 <- ggplot(sub3, aes(Length, Hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Length vs Trans HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

sub5<- rbind(sub2,sub3)
plot_top_rpkm <- ggplot(sub5, aes(RPKM, fill=Mobility)) +
        geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
        ggtitle("Non-Mobile vs Mobile Abundance vs Trans HiC")+
        scale_fill_manual(values = c("orange", "purple")) +
        theme(legend.position = "none")

plot_top_length <- ggplot(sub5, aes(Length, fill=Mobility)) +
  geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
        ggtitle("Non-Mobile vs Mobile Length vs Trans HiC")+
        scale_fill_manual(values = c("orange", "purple")) +
        theme(legend.position = "none")

dense8 <- ggplot(sub2, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense9 <- ggplot(sub3, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

pdf("Densities_scaffolds_all_trans.pdf",height=12,width=10)
grid.arrange(dense2, dense5,dense3, dense6, plot_top_rpkm, plot_top_length, nrow = 3, ncol=2)
dev.off()

pdf("Length_vs_Abundance_scaffolds_all_trans.pdf",height=12,width=5)
grid.arrange(dense8,dense9, nrow=3,ncol=1)
dev.off()
##################################################################################
#cis stuff
sub2 <- subset(subdata, Mobility==0)
dense2 <- ggplot(sub2, aes(RPKM, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Abundance vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

sub3 <- subset(subdata, Mobility==1)
dense3 <- ggplot(sub3, aes(RPKM, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Abundance vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

#density plot of each

dense5 <- ggplot(sub2, aes(Length, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Length vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense6 <- ggplot(sub3, aes(Length, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Length vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

sub5<- rbind(sub2,sub3)
plot_top_rpkm <- ggplot(sub5, aes(RPKM, fill=Mobility)) +
        geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
        ggtitle("Non-Mobile vs Mobile Abundance vs Cis HiC")+
        scale_fill_manual(values = c("orange", "purple")) +
        theme(legend.position = "none")

plot_top_length <- ggplot(sub5, aes(Length, fill=Mobility)) +
  geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
        ggtitle("Non-Mobile vs Mobile Length vs Cis HiC")+
        scale_fill_manual(values = c("orange", "purple")) +
        theme(legend.position = "none")

dense8 <- ggplot(sub2, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense9 <- ggplot(sub3, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

pdf("Densities_scaffolds_all_cis.pdf",height=12,width=10)
grid.arrange(dense2, dense5,dense3, dense6, plot_top_rpkm, plot_top_length, nrow = 3, ncol=2)
dev.off()

pdf("Length_vs_Abundance_scaffolds_all_cis.pdf",height=12,width=5)
grid.arrange(dense8,dense9, nrow=3,ncol=1)
dev.off()
################################################################cis only scaffolds
sub2 <- subset(subdata, Mobility==0 & Cis_hic_hits>0)
dense2 <- ggplot(sub2, aes(RPKM, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Abundance vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

sub3 <- subset(subdata, Mobility==1 & Cis_hic_hits>0)
dense3 <- ggplot(sub3, aes(RPKM, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Abundance vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")


#density plot of each

dense5 <- ggplot(sub2, aes(Length, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Length vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense6 <- ggplot(sub3, aes(Length, Cis_hic_hits))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Length vs Cis HiC")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")


sub5<- rbind(sub2,sub3)
plot_top_rpkm <- ggplot(sub5, aes(RPKM, fill=Mobility)) +
        geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
        ggtitle("Non-Mobile vs Mobile Abundance vs Cis HiC")+
        scale_fill_manual(values = c("orange", "purple")) +
        theme(legend.position = "none")

plot_top_length <- ggplot(sub5, aes(Length, fill=Mobility)) +
  geom_density(alpha=.5) +
        scale_x_log10(labels=comma) +
        ggtitle("Non-Mobile vs Mobile Length vs Cis HiC")+
        scale_fill_manual(values = c("orange", "purple")) +
        theme(legend.position = "none")

dense8 <- ggplot(sub2, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Non-Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("orange",50),high="orange")

dense9 <- ggplot(sub3, aes(Length, RPKM))+
        geom_bin2d()+
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        ggtitle("Mobile Length vs Abundance")+
        scale_fill_continuous(low=makeTransparent("purple",50),high="purple")

pdf("Densities_scaffolds_withcis.pdf",height=12,width=10)
grid.arrange(dense2, dense5,dense3, dense6, plot_top_rpkm, plot_top_length, nrow = 3, ncol=2)
dev.off()

pdf("Length_vs_Abundance_scaffolds_withcis.pdf",height=12,width=5)
grid.arrange(dense8,dense9, nrow=3,ncol=1)
dev.off()

#marginal density of y - plot on the right
#plot_right <- ggplot(subdata, aes(Hic_hits, fill=Mobility)) + 
#  geom_density(alpha=.5) + 
#  coord_flip() + 
#	scale_x_log10() +
#  scale_fill_manual(values = c("orange", "purple")) + 
#  theme(legend.position = "none") 

sub <- subset(subdata, Hic_hits>0)
cor.test(sub$RPKM, sub$Hic_hits, method="spearman")
cor.test(sub$Length, sub$Hic_hits, method="spearman")
sub <- subset(subdata, Hic_hits>0 & Mobility==0)
cor.test(sub$RPKM, sub$Hic_hits, method="spearman")
cor.test(sub$Length, sub$Hic_hits, method="spearman")
sub <- subset(subdata, Hic_hits>0 & Mobility ==0)
cor.test(sub$RPKM, sub$Hic_hits, method="spearman")
cor.test(sub$Length, sub$Hic_hits, method="spearman")
sub <- subset(subdata, Hic_hits>0 & Mobility == 1)
cor.test(sub$RPKM, sub$Hic_hits, method="spearman")
cor.test(sub$Length, sub$Hic_hits, method="spearman")

sub <- subset(subdata, Cis_hic_hits>0)
cor.test(sub$RPKM, sub$Cis_hic_hits, method="spearman")
cor.test(sub$Length, sub$Cis_hic_hits, method="spearman")
sub <- subset(subdata, Cis_hic_hits>0 & Mobility==0)
cor.test(sub$RPKM, sub$Cis_hic_hits, method="spearman")
cor.test(sub$Length, sub$Cis_hic_hits, method="spearman")
sub <- subset(subdata, Cis_hic_hits>0 & Mobility ==0)
cor.test(sub$RPKM, sub$Cis_hic_hits, method="spearman")
cor.test(sub$Length, sub$Cis_hic_hits, method="spearman")
sub <- subset(subdata, Cis_hic_hits>0 & Mobility == 1)
cor.test(sub$RPKM, sub$Cis_hic_hits, method="spearman")
cor.test(sub$Length, sub$Cis_hic_hits, method="spearman")

sub <- subset(subdata, Hic_hits>0)
cor.test(sub$RPKM, sub$Hic_hits)
cor.test(sub$Length, sub$Hic_hits)
sub <- subset(subdata, Hic_hits>0 & Mobility==0)
cor.test(sub$RPKM, sub$Hic_hits)
cor.test(sub$Length, sub$Hic_hits)
sub <- subset(subdata, Hic_hits>0 & Mobility ==0)
cor.test(sub$RPKM, sub$Hic_hits)
cor.test(sub$Length, sub$Hic_hits)
sub <- subset(subdata, Hic_hits>0 & Mobility == 1)
cor.test(sub$RPKM, sub$Hic_hits)
cor.test(sub$Length, sub$Hic_hits)

sub <- subset(subdata, Cis_hic_hits>0)
cor.test(sub$RPKM, sub$Cis_hic_hits)
cor.test(sub$Length, sub$Cis_hic_hits)
sub <- subset(subdata, Cis_hic_hits>0 & Mobility==0)
cor.test(sub$RPKM, sub$Cis_hic_hits)
cor.test(sub$Length, sub$Cis_hic_hits)
sub <- subset(subdata, Cis_hic_hits>0 & Mobility ==0)
cor.test(sub$RPKM, sub$Cis_hic_hits)
cor.test(sub$Length, sub$Cis_hic_hits)
sub <- subset(subdata, Cis_hic_hits>0 & Mobility == 1)
cor.test(sub$RPKM, sub$Cis_hic_hits)
cor.test(sub$Length, sub$Cis_hic_hits)
#################################ok now just do it without subsetting
cor.test(subdata$RPKM, subdata$Hic_hits, method="spearman")
cor.test(subdata$Length, subdata$Hic_hits, method="spearman")
sub <- subset(subdata, Mobility==0)
cor.test(sub$RPKM, sub$Hic_hits, method="spearman")
cor.test(sub$Length, sub$Hic_hits, method="spearman")
sub <- subset(subdata, Mobility ==0)
cor.test(sub$RPKM, sub$Hic_hits, method="spearman")
cor.test(sub$Length, sub$Hic_hits, method="spearman")
sub <- subset(subdata, Mobility == 1)
cor.test(sub$RPKM, sub$Hic_hits, method="spearman")
cor.test(sub$Length, sub$Hic_hits, method="spearman")

cor.test(subdata$RPKM, subdata$Cis_hic_hits, method="spearman")
cor.test(subdata$Length, subdata$Cis_hic_hits, method="spearman")
sub <- subset(subdata, Mobility==0)
cor.test(sub$RPKM, sub$Cis_hic_hits, method="spearman")
cor.test(sub$Length, sub$Cis_hic_hits, method="spearman")
sub <- subset(subdata, Mobility ==0)
cor.test(sub$RPKM, sub$Cis_hic_hits, method="spearman")
cor.test(sub$Length, sub$Cis_hic_hits, method="spearman")
sub <- subset(subdata, Mobility == 1)
cor.test(sub$RPKM, sub$Cis_hic_hits, method="spearman")
cor.test(sub$Length, sub$Cis_hic_hits, method="spearman")

#pearson
cor.test(subdata$RPKM, subdata$Hic_hits)
cor.test(subdata$Length, subdata$Hic_hits)
sub <- subset(subdata, Mobility==0)
cor.test(sub$RPKM, sub$Hic_hits)
cor.test(sub$Length, sub$Hic_hits)
sub <- subset(subdata, Mobility ==0)
cor.test(sub$RPKM, sub$Hic_hits)
cor.test(sub$Length, sub$Hic_hits)
sub <- subset(subdata, Mobility == 1)
cor.test(sub$RPKM, sub$Hic_hits)
cor.test(sub$Length, sub$Hic_hits)

cor.test(subdata$RPKM, subdata$Cis_hic_hits)
cor.test(subdata$Length, subdata$Cis_hic_hits)
sub <- subset(subdata, Mobility==0)
cor.test(sub$RPKM, sub$Cis_hic_hits)
cor.test(sub$Length, sub$Cis_hic_hits)
sub <- subset(subdata, Mobility ==0)
cor.test(sub$RPKM, sub$Cis_hic_hits)
cor.test(sub$Length, sub$Cis_hic_hits)
sub <- subset(subdata, Mobility == 1)
cor.test(sub$RPKM, sub$Cis_hic_hits)
cor.test(sub$Length, sub$Cis_hic_hits)


#######################################
library(scales)
library(plyr)
indata$hashic<-indata$Hic_hits>0
indata$hascishic <-indata$Cis_hic_hits>0
indata$count<- rep(1,length(indata$hashic))

mu_len <- ddply(indata, "hashic", summarise, grp.mean=median(Length))
mu_abund<-ddply(indata, "hashic", summarise, grp.mean=median(RPKM))

mu_cis_len <- ddply(indata, "hascishic", summarise, grp.mean=median(Length))
mu_cis_abund<-ddply(indata, "hascishic", summarise, grp.mean=median(RPKM))

pdf("Contigs_with_without_hiclinks.pdf")
ggplot(data=indata,aes(hashic))+
geom_bar()+
scale_y_continuous(labels = comma)
dev.off()

pdf("Contigs_with_without_cishiclinks.pdf")
ggplot(data=indata,aes(hascishic))+
geom_bar()+
scale_y_continuous(labels = comma)
dev.off()

pdf("Contigs_with_without_hiclinks_vs_length.pdf")
ggplot(indata, aes(x=Length,fill=hashic, color=hashic)) +
	geom_histogram( position="dodge")+
	scale_x_log10(labels=comma)+
	scale_fill_manual(name="Contig Has Trans Reads",values=c("blue","red"))+
	scale_color_manual(name="Median Length", values=c("blue","red"))+
	geom_vline(data=mu_len, aes(xintercept=grp.mean, color=hashic),linetype="dashed")
dev.off()

pdf("Contigs_with_without_cishiclinks_vs_length.pdf")
ggplot(indata, aes(x=Length,fill=hascishic, color=hascishic)) +
        geom_histogram(position="dodge")+
        scale_x_log10(labels=comma)+
	ylab("Number of Contigs")+
	scale_fill_manual(name="Contig Has Cis Reads",values=c("blue","red"))+
        scale_color_manual(name="Median Length", values=c("blue","red"))+
        geom_vline(data=mu_cis_len, aes(xintercept=grp.mean, color=hascishic),linetype="dashed")
dev.off()

pdf("Contigs_with_without_hiclinks_vs_abundance.pdf")
ggplot(indata, aes(x=RPKM,fill=hashic, color=hashic)) +
	geom_histogram(position="dodge")+
	scale_x_log10(labels=comma)+
	ylab("Number of Contigs")+
        scale_fill_manual(name="Contig Has Trans Reads",values=c("blue","red"))+
        scale_color_manual(name="Median Abundance", values=c("blue","red"))+	
	geom_vline(data=mu_abund, aes(xintercept=grp.mean, color=hashic),linetype="dashed")
dev.off()


pdf("Contigs_with_without_cishiclinks_vs_abundance.pdf")
ggplot(indata, aes(x=RPKM,fill=hascishic, color=hascishic)) +
        geom_histogram(position="dodge")+
        scale_x_log10(labels=comma)+
        ylab("Number of Contigs")+
        scale_fill_manual(name="Contig Has Cis Reads",values=c("blue","red"))+
        scale_color_manual(name="Median Abundance", values=c("blue","red"))+      
        geom_vline(data=mu_cis_abund, aes(xintercept=grp.mean, color=hascishic),linetype="dashed")
dev.off()
############


########################################################################
#looking at abundance
########################################################################

pdf("Density_length_vs_hictranshits_bin.pdf")
ggplot(indata, aes(Length, Hic_hits)) +
        geom_bin2d() +
        scale_x_continuous(labels=comma)+
        scale_y_continuous(labels=comma)+
        scale_fill_viridis() 
dev.off()
###############################################################################################
sub <- subset(indata, Length>0 & Hic_hits>0)
sub <- subset(indata, Length>0 & Hic_hits>0 & Mobility==1)
pdf("Density_length_vs_hictranshits_bin_subset_mobile.pdf")
ggplot(sub, aes(Length, Hic_hits)) +
        geom_bin2d(bins=50) +
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        scale_fill_viridis()
dev.off()

sub <- subset(indata, Length>0 & Hic_hits>0 & Mobility==0)
pdf("Density_length_vs_hictranshits_bin_subset_nonmobile.pdf")
ggplot(sub, aes(Length, Hic_hits)) +
        geom_bin2d(bins=50) +
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        scale_fill_viridis()
dev.off()

sub <- subset(indata, Length>0 & Hic_hits>0 & Mobility==0)
pdf("Density_length_vs_hictranshits_bin_subset_nonmobile.pdf")
ggplot(sub, aes(Length, Hic_hits)) +
        geom_bin2d(bins=50) +
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)+
        scale_fill_viridis()
dev.off()







##########hic cis hits
###############################################
pdf("Density_length_vs_hiccishits_bin.pdf")
ggplot(indata, aes(Length, Cis_hic_hits)) +
        geom_bin2d() +
        scale_x_continuous(labels=comma)+
        scale_y_continuous(labels=comma)+
        scale_fill_viridis()
dev.off()

pdf("Density_length_vs_hiccishits_bin_loglog.pdf")
ggplot(indata, aes(Length, Cis_hic_hits)) +
        geom_bin2d() +
        #scale_fill_viridis() +
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)
dev.off()

pdf("Density_abundance_vs_hiccishits_bin.pdf")
ggplot(indata, aes(RPKM, Cis_hic_hits)) +
  geom_bin2d() +
        scale_x_continuous(labels=comma)+
        scale_y_continuous(labels=comma)
  scale_fill_viridis()
dev.off()

pdf("Density_abundance_vs_hiccishits_bin_loglog.pdf")
ggplot(indata, aes(RPKM, Cis_hic_hits)) +
  geom_bin2d() +
  scale_fill_viridis() +
        scale_x_log10(labels=comma)+
        scale_y_log10(labels=comma)
dev.off()

