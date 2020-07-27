#map the names back #this is specific to 95%
arg_name_dict = {}
prot_arg_name = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95_names.txt'
with open(prot_arg_name) as f:
	for line in f:
		r = line.split('\t')[0].strip()
		if r != 'ORF_ID':
			arg_name_dict[r] = (line.strip().split('\t')[2], line.strip().split('\t')[3])


#make a shortened list of argnames
argnames = list(set(arg_name_dict.values()))
argnames.sort()

inhandle = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes3/args_' + '95' + '_nr.fna.clstr'
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

listofcoinames = []
incoifile = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/COI_sort.txt'
with open(incoifile) as incoi:
	for line in incoi:
		listofcoinames.append(line.split(':')[0])
		listofcoinames.append(line.split(':')[1].strip())

coinames = list(set(listofcoinames))		
coinames.remove('NA')

outhandle = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/COI_clusters.txt"
with open(outhandle,'w') as outfile:
	for key in clusterprotdict.keys():
		genes = clusterprotdict[key]
		resname = arg_name_dict[clusternum_map[key]][0]
		cardname=arg_name_dict[clusternum_map[key]][1]
		if ((resname in coinames) or (cardname in coinames)):
			outfile.write(key + '\t' + resname + '\t' + cardname + '\n')
