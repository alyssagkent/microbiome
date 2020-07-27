library (RColorBrewer)
library(gplots)
library(ggplot2)
library(tsne)
library(vegan)
library(reshape)
###ggplot theme no background
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
  theme(panel.grid = element_blank())}
theme_set(theme_nogrid())


#########set variable things

#/programs/MetaPhlAn-2.0/utils/merge_metaphlan_tables.py *_profile.txt > CDC+Healthy+Press_metaphlan.txt
folder='/workdir/users/agk85/CDC2/metaphlan/cdc'
infile='/workdir/users/agk85/CDC2/metaphlan/cdc/mgm_hic_metaphlan.txt'
pairfile = '/workdir/users/agk85/CDC2/metaphlan/cdc/mgm_hic_pairs_cdc+healthy+press.csv'
generalname = 'CDC+Healthy+Press'

#make sure it doesn't ahve profile stuff
command = paste("sed -i 's/_profile//g' ",infile, sep="") 
system(command)

keepers <- c("B370.1","B370.1hic","B370.2","B370.2hic","B370.3","B370.3hic",
	"B370.4","B370.4hic","B370.5","B370.5hic",
        "B357.1","B357.1hic","B357.2","B357.2hic","B357.3","B357.3hic",
        "B357.4","B357.4hic","B357.5","B357.5hic","B314.1","B314.1hic",
        "B314.2","B314.2hic","B314.3","B314.3hic","B314.4","B314.4hic",
        "B316.1","B316.1hic","B316.2","B316.2hic","B316.3","B316.3hic",
        "B316.4","B316.4hic","B316.5","B316.5hic","B316.6","B316.6hic",
        "B316.7","B316.7hic","B320.1","B320.1hic","B320.2","B320.2hic",
        "B320.3","B320.3hic","B320.5","B320.5hic",
        "B331.1","B331.1hic","B331.2","B331.2hic","B331.3","B331.3hic",
        "B331.4","B331.4hic","B335.1","B335.1hic","B335.2","B335.2hic",
        "B335.3","B335.3hic","B357.6","B357.6hic",
	"ProxiMeta.1","MluCI.1hic","Sau3aI.1hic",
	"US3.8","US3.10","US3.12","US3.14","US3.16",
	"US3.8hic","US3.10hic","US3.12hic","US3.14hic","US3.16hic",
	"US8.1","US8.2","US8.3","US8.4","US8.5",
	"US8.1hic","US8.2hic","US8.3hic","US8.4hic","US8.5hic")
#########
setwd(folder)
metaphlan1 = read.csv(infile, sep="\t",header=T,stringsAsFactors = FALSE)
#remove the weird id line
metaphlan <- metaphlan1[!grepl("#SampleID", metaphlan1$ID),]
#levels k=1,p=2,c=3,o=4,f=5,g=6,s=7,t=8
#set your level:
taxnames = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
levs = c(2,3,4,5,6,7)
for (l in levs){
taxalev=taxnames[l]
levels <-as.vector(sapply(strsplit(as.character(metaphlan$ID),"\\|"), function(x) length(x)))
metaphlan$levels <- levels
metaphlan_sub <- subset(metaphlan, levels == l)

#this gets a single name from the ID
taxa <-sapply(strsplit(as.character(metaphlan_sub$ID),"__"), `[`, length(strsplit(as.character(metaphlan_sub$ID),"__")[[l+1]]))
rownames(metaphlan_sub)<-taxa
#this removes a few columns before making it a matrix
df<-metaphlan_sub[,!colnames(metaphlan_sub) %in% c("levels","ID","Blank_profile")]
Patient<-as.vector(sapply(strsplit(as.character(colnames(df)),"\\."), function(x) x[1]))
charlist<- Patient # colnames(df)
colors = rainbow(length(unique(charlist)))
names(colors) = unique(charlist)
ecb = function(x,y){ plot(x,t='n'); text(x,labels=colnames(df),cex=.5, col=colors[charlist]) }

#subset to the libraries of interest
d <- data.frame(data.matrix(df[,colnames(df) %in% keepers]))

pdf(paste(generalname , "_hclust_vegdist_",taxalev,".pdf",sep=""))
 plot(hclust(vegdist(t(d))), cex=.35)
 dev.off()

write.table(as.matrix(vegdist(t(d))), paste(generalname,"_", taxalev, "_vegdist.txt",sep=""), row.names=T, sep = '\t')

mgm_hic<-read.csv(pairfile,header=T)
mgm_hic_cor <- data.frame(Name=character(),
                 r=numeric(),p=numeric()) 
all_mgm_hic <- data.frame(n=character(),
                 m=numeric(),h=numeric()) 
all_mgm_distlist <- data.frame(m=numeric(),h=numeric())
 
for (i in 1:(nrow(mgm_hic))){
mgm = as.character(mgm_hic[i,1])
hic = as.character(mgm_hic[i,2])

m = d[[mgm]]
h = d[[hic]]
mhdist = vegdist(t(c(m,h)))
r = cor.test(m,h)$estimate
p = cor.test(m,h)$p.value
v = data.frame(hic, r, p)
mgm_hic_cor<-rbind(mgm_hic_cor, v)
Patient<-strsplit(mgm,"\\.")[[1]][1]
sp <- data.frame(m, h)
rownames(sp)<-rownames(d)
sp_clean <- sp[rowSums(sp)>0,]
sp_clean$r <- rownames(sp_clean)
sp_clean$sample <- hic
sp_clean$patient <- strsplit(hic, '[.]')[[1]][1]
sp_clean$names <- paste(rownames(sp_clean), rep(mgm, nrow(sp_clean)), sep="_")
all_mgm_hic <- rbind(all_mgm_hic, sp_clean)	
}
all_mgm_hic$diff <- (all_mgm_hic$h - all_mgm_hic$m)
all_mgm_hic_pre <- subset(all_mgm_hic, patient %in% c("B314","B316","B320","B331","B335","B357","B370"))
all_mgm_hic_pre$prepost <- rep("CDC", nrow(all_mgm_hic_pre))

all_mgm_hic_post <- subset(all_mgm_hic, patient %in% c("Sau3aI","MluCI","ProxiMeta"))
all_mgm_hic_post$prepost <- rep("Press", nrow(all_mgm_hic_post))

all_mgm_hic_healthy <- subset(all_mgm_hic, patient %in% c("US3","US8"))
all_mgm_hic_healthy$prepost <-rep("Healthy",nrow(all_mgm_hic_healthy))

all_mgm_hic_both<- rbind(all_mgm_hic_pre,all_mgm_hic_post,all_mgm_hic_healthy)

#Differences
pdf(paste(generalname, "_metagenome_differences_",taxalev,".pdf",sep=""), height=8, width = 5*l,useDingbats=FALSE)
g<-ggplot(data=all_mgm_hic_both, aes(r,diff))+
	geom_boxplot(outlier.shape = NA)+
	geom_point(aes(col=patient, shape=prepost),position=position_jitter())+
	ylab("Relative Abundance Difference (HiC - Metagenome)")+ xlab(taxalev)+
	ylim(-100,100)+
	scale_shape_manual("Restriction Enzyme Round",labels=c("CDC"="CDC", "Press"="ProxiMeta", "Healthy"="Healthy"), values=c(19,17,15))+
	scale_color_manual("Patient", values=c("blue","mediumblue","midnightblue","skyblue3","steelblue3","turquoise2","navy","violetred2","red","orange","sienna3","sienna4","sandybrown","tomato"))+
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(g) 
dev.off()




a = subset(all_mgm_hic_both, patient %in% c("Sau3aI","MluCI","ProxiMeta"))
pdf(paste("Proximeta", generalname, "_metagenome_differences_",taxalev,".pdf",sep=""), height=8, width = 5*l,useDingbats=FALSE)
g<-ggplot(data=a, aes(r,diff))+
        geom_boxplot(outlier.shape = NA)+
        geom_point(aes(col=patient, shape=prepost),position=position_jitter())+
        ylab("Relative Abundance Difference (HiC - Metagenome)")+ xlab(taxalev)+
        ylim(-100,100)+
        scale_shape_manual("Restriction Enzyme Round",labels=c("CDC"="CDC", "Press"="ProxiMeta", "Healthy"="Healthy"), values=c(19,17,15))+
        scale_color_manual("Patient", values=c("blue","mediumblue","midnightblue","skyblue3","steelblue3","turquoise2","navy","violetred2","red","orange","sienna3","sienna4","sandybrown","tomato"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(g)
dev.off() 
####stacked bar chart hic and metagenomes
colourCount <- length(unique(all_mgm_hic$r))
all_mgm_hic$diff<-NULL
all_mgm_hic$diffnorm<-NULL
all_mgm_hic_melt <- melt(all_mgm_hic)
pdf(paste(generalname, "_stacked_",taxalev,".pdf",sep=""), height = 5*l, width=30,useDingbats=FALSE)
g<-ggplot(data=all_mgm_hic_melt, aes(variable,value, fill=r))+
	geom_col(position="stack")+
	xlab("Metagenome and Hi-C Samples")+
	ylab(paste(taxalev, " Relative Abundance", sep=""))+
	labs(fill='Taxa')+
	scale_fill_manual(taxalev, values = colorRampPalette(brewer.pal(8,"Accent"))(colourCount))+
	theme(legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1))+
	facet_grid(.~sample)
print(g)
dev.off()
}

