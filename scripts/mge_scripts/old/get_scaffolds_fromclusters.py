
scfdict = {}
with open('mge_95_nr.fna.clstr.tbl.annot') as reffile:
	for line in reffile:
		clusterid = line.split('\t')[0]
		pfam = line.split('\t')[1]
		rep = line.split('\t')[3]
		genes = line.split('\t')[4]
		scfdict[clusterid] = (rep, genes, pfam)

inhandle = 'cluster_phage_shared.txt'
with open(inhandle) as infile:
	for line in infile:
		clusterid = line.split('_')[1]
		rep, genes, pfam = scfdict[clusterid]
		print(clusterid)
		print(genes)
		print(pfam)


