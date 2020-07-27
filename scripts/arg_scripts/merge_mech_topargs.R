df1 <- read.csv('/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech.txt',header=T, sep="\t")
df2 <- read.csv('/workdir/users/agk85/CDC2/args/patric_comparisons/all_clusters_topargs_mergednames.txt',header=T,sep="\t")

dfmerge <- merge(df1, df2, by=c("Cluster"), all.x=T)
write.table(dfmerge,"/workdir/users/agk85/CDC2/args/toparg_mechs.txt",quote=F,row.names=F,sep="\t")
