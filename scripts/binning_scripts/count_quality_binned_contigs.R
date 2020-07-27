
bins = read.table('~/agk/CDC2/das/all_DASTool_scaffolds2bin.txt')
colnames(bins) <- c("Contigs","Bins")
bintables=read.table('~/agk/CDC2/das/all_bintables.txt')
colnames(bintables) <- c("Bins","Patient","Sample","Quality","Taxonomy")

mergedf = merge(bins, bintables, by=c("Bins"))
