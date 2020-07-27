library(ggplot2)
library(scales)
library(Hmisc)
#Rscript
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())
setwd("/workdir/users/agk85/CDC2/args/min_dist")

args <- commandArgs(trailingOnly = TRUE)
minreads=args[1]
genetype=args[2]



infile <- read.csv(paste(genetype, "_dist_to_nearest_hicread_min_",minreads,".txt",sep=""),sep="\t", header=T)
infile$Ratio <- infile$Minimum_distance/infile$Scaffold_length
ggplot(data=infile,aes(fill = factor(Scf_taxonomy)))+
geom_histogram(aes(Ratio),position="stack")

ggplot(data=infile)+
geom_density2d(aes(Scaffold_length, Minimum_distance))+
scale_y_log10(labels=comma)+
scale_x_log10(labels=comma)


pdf(paste("Min_dist_vs_scaffold_density_",genetype,"_",minreads,".pdf",sep=""))
ggplot(infile, aes(x = Scaffold_length,y =Minimum_distance))+
stat_binhex()+
scale_y_continuous(labels=comma)+
scale_x_continuous(labels=comma)
dev.off()

pdf(paste("Min_dist_vs_scaffold_density_loglog_",genetype,"_",minreads,".pdf",sep=""))
ggplot(infile, aes(x = Scaffold_length,y =Minimum_distance))+
stat_binhex()+
scale_y_log10(labels=comma,limits=c(1,1000000))+
scale_x_log10(labels=comma,limits=c(1,1000000))
dev.off()

sub <- subset(infile, Scf_taxonomy != 1)

#which ones equal 0?? that means they have a read hitting INSIDDE the  ARG
subzero <- subset(sub, Minimum_distance==0)

#make the histogram, but split into groups by the Scaffold_length
#groups 0-10,000, 10,000 to 50,000; 50,000, 


pdf(paste("Min_dist_gene_nogroups_ratio_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub)+
geom_histogram(aes(Ratio), bins=50)+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()

#split into ratio groups
s1 <- subset(sub, Ratio <=.5)
s2 <- subset(sub, Ratio >0.5)






sub$groups <- cut2(sub$Scaffold_length, c(10000,50000,100000,500000,1000000)) 

pdf(paste("Min_dist_gene_groupcolor_distance_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub, aes(fill=groups))+
geom_histogram(aes(Minimum_distance), position="stack",binwidth=10000)+
scale_x_continuous(labels=comma)+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()


pdf(paste("Min_dist_gene_groupcolor_distance_inset_0_100000_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub, aes(fill=groups))+
geom_histogram(aes(Minimum_distance), position="stack",bins=50)+
scale_x_continuous(labels=comma, limits=c(0,100000))+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()


pdf(paste("Min_dist_gene_groupcolor_distance_inset_0_50000_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub, aes(fill=groups))+
geom_histogram(aes(Minimum_distance), position="stack",bins=50)+
scale_x_continuous(labels=comma, limits=c(0,50000))+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()

pdf(paste("Min_dist_gene_groupcolor_distance_inset_0_20000_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub, aes(fill=groups))+
geom_histogram(aes(Minimum_distance), position="stack",bins=50)+
scale_x_continuous(labels=comma, limits=c(0,20000))+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()





#####################################OLD STUFF

pdf(paste("OLD_Min_dist_gene_groupcolor_ratio_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub, aes(fill=groups))+
geom_histogram(aes(Ratio), position="stack",bins=50)+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()

pdf(paste("OLD_Min_dist_gene_groupcolor_distance_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub, aes(fill=groups))+
geom_histogram(aes(Minimum_distance), position="stack",bins=50)+
scale_x_continuous(labels=comma)+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()

pdf(paste("OLD_Min_dist_gene_groupcolor_distance_loglog_",genetype,"_",minreads,".pdf",sep=""))
ggplot(data=sub, aes(fill=groups))+
geom_histogram(aes(Minimum_distance), position="stack",bins=50)+
scale_x_log10(labels=comma)+
ylab("Number of Gene-Taxonomy Combinations")
dev.off()

pdf(paste("OLD_Min_dist_gene_groups_distance_",genetype,"_",minreads,".pdf",sep=""),height=7, width=20)
ggplot(data=sub)+
geom_histogram(aes(Minimum_distance))+
xlab("Min. Distance to Nearest Read")+
facet_grid(.~groups, scales="free")+
scale_x_continuous(labels=comma)
dev.off()

#normalize it
pdf(paste("OLD_Min_dist_gene_groups_ratio_",genetype,"_",minreads,".pdf",sep=""), height=7, width = 20)
ggplot(data=sub)+
geom_histogram(aes(Ratio))+
xlab("Ratio (Min. Distance to Nearest Read/Scaffold Length)")+
facet_grid(.~groups, scales="free")+
scale_x_continuous(labels=comma)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


#so if you want to aggregate down to the min distance dist of each patient-ARG cluster, then you aggregate
sub.agg <- aggregate(Minimum_distance ~ Cluster+Patient + Taxonomy, data=sub, min)
sub.trim <- subset(sub, select=c("Cluster","Patient","Taxonomy","Minimum_distance","Scaffold_length"))
#but this doesn't get you the scaffold length so add it back in
#but that means you sometimes have multiple scaffolds to one combination of min-dist
sub.agg.add <- merge(sub.agg, sub.trim, by=c("Cluster","Patient","Taxonomy","Minimum_distance"))
sub.agg.add$groups <- cut2(sub.agg.add$Scaffold_length, c(10000,50000,100000,500000,1000000)) 
pdf(paste("OLD_Min_dist_argpat_groups_distance_",genetype,"_",minreads,".pdf",sep=""),height=7, width=20)
ggplot(data=sub.agg.add)+
geom_histogram(aes(Minimum_distance))+
xlab("Min. Distance to Nearest Read")+
facet_grid(.~groups, scales="free")+
scale_x_continuous(labels=comma)
dev.off()


#so if you want to aggregate down to the min distance dist of each patient-ARG cluster, then you aggregate
sub.agg <- aggregate(Ratio ~ Cluster+Patient + Taxonomy, data=sub, min)
#normalize it
pdf(paste("OLD_Min_dist_argpat_ratio_",genetype,"_",minreads,".pdf",sep=""), height=7, width = 20)
ggplot(data=sub.agg)+
geom_histogram(aes(Ratio))+
xlab("Ratio (Min. Distance to Nearest Read/Scaffold Length)")+
scale_x_continuous(labels=comma)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#########################
#ok now do it on the infile data so with scaffolds that have taxonomies
#so if you want to aggregate down to the min distance dist of each patient-ARG cluster, then you aggregate
sub.agg <- aggregate(Minimum_distance ~ Cluster+Patient + Taxonomy, data=infile, min)
sub.trim <- subset(infile, select=c("Cluster","Patient","Taxonomy","Minimum_distance","Scaffold_length"))
#but this doesn't get you the scaffold length so add it back in
#but that means you sometimes have multiple scaffolds to one combination of min-dist
sub.agg.add <- merge(sub.agg, sub.trim, by=c("Cluster","Patient","Taxonomy","Minimum_distance"))
sub.agg.add$groups <- cut2(sub.agg.add$Scaffold_length, c(10000,50000,100000,500000,1000000)) 
pdf(paste("OLD_Min_dist_argpat_groups_distance_withscftaxa_",genetype,"_",minreads,".pdf",sep=""),height=7, width=20)
ggplot(data=sub.agg.add)+
geom_histogram(aes(Minimum_distance))+
xlab("Min. Distance to Nearest Read")+
facet_grid(.~groups, scales="free")+
scale_x_continuous(labels=comma)
dev.off()


#so if you want to aggregate down to the min distance dist of each patient-ARG cluster, then you aggregate
sub.agg <- aggregate(Ratio ~ Cluster+Patient + Taxonomy, data=infile, min)
#normalize it
pdf(paste("OLD_Min_dist_argpat_ratio_withscftaxa_",genetype,"_",minreads,".pdf",sep=""), height=7, width = 20)
ggplot(data=sub.agg)+
geom_histogram(aes(Ratio))+
xlab("Ratio (Min. Distance to Nearest Read/Scaffold Length)")+
scale_x_continuous(labels=comma)+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



