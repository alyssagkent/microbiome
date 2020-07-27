# blast_link_scaffold_reference.R
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)

setwd("~/agk/synthetic/reference_genomes/")

blastout <- read.csv("Synthetic-50_scaffold.out", header=F, sep="\t")
header=c("Query","Subject","PID","Length")
blast <- blastout[,1:4]
colnames(blast)<-header

#aggregate on the first two columns summing on the fourth
#create average percent identity and sum the lengths across

blast$Good <- blast$PID*blast$Length
by_query_subject<-group_by(blast, Query,Subject)

blastagg <- summarise(by_query_subject,n = n(),good=sum(Good, na.rm=TRUE),length = sum(Length, na.rm = TRUE))

blastagg$avg_pid <- blastagg$good/blastagg$length

#max of the blastagg based on only the Query --- aka which genome matches more
a= data.frame(blastagg %>% group_by(Query) %>% top_n(1, length))
by_query <-group_by(a, Query)

#from these take the top query based on higher average percent identity
b= data.frame(by_query %>% group_by(Query) %>% top_n(1, avg_pid))

#figure out which ones have more than 2 still yet
c <- data.frame(summarise(by_query,n = n(),ap=mean(avg_pid, na.rm=TRUE),l = sum(length, na.rm = TRUE)))

#subset to those multi-mapping contigs
multi <- subset(c, n>1)
#make a list of the multi_mapping file
multi_mapping <- multi$Query

#remove the multi_mapping files
b_nomulti <- subset(b, !(Query %in% multi_mapping))

#write the association file
write.table(b_nomulti, "Synthetic-50_scaffold_bestguess.txt",sep="\t",quote=FALSE,row.names=FALSE)



