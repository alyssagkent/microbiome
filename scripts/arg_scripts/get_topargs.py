#output clusters that have top genes
#load the card aro ids that are top genes
#two files to get beforehand
# cat */*resfams.tbl.txt cat {output.restable}  | grep -v '^#' | awk '{{print $1,$3,$4}}' | sed 's/ /\t/g' > all_resfams_rf.txt

# cat */*resfams.tbl.txt cat {output.restable}  | grep -v '^#' | awk '{{print $1,$3,$4}}' | sed 's/ /\t/g' > all_resfams_rf.txt
# 

inhandle = '/workdir/users/agk85/CDC2/args/args_99_nr.fna.clstr'
protclusterdict = {}
clusterprotdict = {}
clusternum_map = {}
with open(inhandle) as infile:
	for line in infile:
		if line[0] == '>':
			cluster = line.strip().split(' ')[1]
		else:
			prot =  line.split('>')[1].split('...')[0]
			protclusterdict[prot] = cluster
			try:
				clusterprotdict[cluster].append(prot)
			except:
				clusterprotdict[cluster] = [prot]
			if '*' in line:	
				clusternum_map[cluster] = prot

topargs = []
with open('/workdir/users/agk85/CDC2/args/patric_comparisons/card_aros_topgenes_qnrb.txt') as cardfile:
	for line in cardfile:
		topargs.append(line.strip())

with open('/workdir/users/agk85/CDC2/args/patric_comparisons/resfams_rfs_topgenes.txt') as resfamsfile:
	for line in resfamsfile:
		topargs.append(line.strip())

toparg_clusters = {}
cardhandle = '/workdir/users/agk85/CDC2/args/patric_comparisons/all_card_aros.txt'
resfamshandle = '/workdir/users/agk85/CDC2/args/patric_comparisons/all_resfams_rfs.txt'


with open(resfamshandle) as ref1:
	for line in ref1:
		gene = line.split('\t')[0]
		argname = line.split('\t')[1]
		rfid = line.strip().split('\t')[2]
		if rfid in topargs:
			cluster = protclusterdict[gene]
			try:
				toparg_clusters[cluster].append(argname) 
			except:
				toparg_clusters[cluster] = [argname] 

with open(cardhandle) as ref2:
	for line in ref2:
		gene = line.split('\t')[0].split(' ')[0]
		argname = line.split('\t')[1]
		aroids= line.strip().split('\t')[2].split(',')
		for aroid in aroids:
			if aroid in topargs:
				cluster = protclusterdict[gene]
				try:
					toparg_clusters[cluster].append(argname)
				except:
					toparg_clusters[cluster] = [argname] 


typedict = {}
with open('/workdir/users/agk85/CDC2/args/patric_comparisons/top_args_categories.csv') as typefile:
	for line in typefile:
		print(line)
		argname,argtype,confidence,specific_name,mechanism,broad_mechanism = line.strip().split(',')
		for arg in argname.split(':'):
			if arg != 'NA':
				typedict[arg]= argtype

print(typedict)

outhandle = '/workdir/users/agk85/CDC2/args/patric_comparisons/all_clusters_topargs_mergednames.txt'
with open(outhandle,'w') as outfile:
	outfile.write('Cluster\tMergedname\tAlltypes\n')
	for cluster in toparg_clusters:
		alltypes = []
		mergename = '|'.join(list(set(toparg_clusters[cluster])))
		for genename in list(set(toparg_clusters[cluster])):
			try:
				alltypes.append(typedict[genename])
			except KeyError:
				print(genename) 
				print("not in top_args_categories.csv")
		outfile.write('{0}\t{1}\t{2}\n'.format(cluster,mergename,','.join(list(set(alltypes)))))
