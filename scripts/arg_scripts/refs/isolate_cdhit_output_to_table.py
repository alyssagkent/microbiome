#make a more parse-able table from cd-hit output
#real specific right now for the mge's
import sys

inhandle = sys.argv[1]
outhandle = sys.argv[2]

protclusterdict = {}
clusterprotdict = {}
clusternum_map = {}
with open(inhandle) as infile:
	for line in infile:
		if line[0] == '>':
			cluster = line.strip().split(' ')[1]
		else:
			prot =  line.split('>')[1].split('...')[0]
			samp = prot.split('_')[0]
			protclusterdict[prot] = cluster
			try:
					clusterprotdict[cluster].append(prot)
			except:
					clusterprotdict[cluster] = [prot]
			if '*' in line:
				clusternum_map[cluster] = prot

isolate_orgs = {}
with open('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/isolate_comparisons/isolate_organisms.txt') as infile:
	for line in infile:
		isolate = line.strip().split('\t')[0]
		org = line.strip().split('\t')[1]
		isolate_orgs[isolate] = org


argdict = {}
prot_arg_name = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95_names.txt'
with open(prot_arg_name) as f:
	for line in f:
		r = line.split('\t')[0].strip()
		if r != 'ORF_ID':
			argdict[r] = line.strip().split('\t')[1]

#make a tab delimited file in this format
name = 'cluster\targname\tstarredid\tids\tnumber_of_proteins\torgs\tdistinct_isolate_species\n'
with open(outhandle, 'w') as outfile:
	for key in clusterprotdict.keys():
		prots = clusterprotdict[key]
		starred = clusternum_map[key]
		try:
			argname = argdict[starred]
		except:
			argname = 'NA'
		orgs = ''
		orglist = []
		for prot in clusterprotdict[key]:
			if 'NODE' in prot:
				isolate = prot.split('_')[0]
				org = isolate_orgs[isolate]
				orgs = orgs + ',' + org
				orglist.append(org)
		num_isolates = str(len(set(orglist)))	
		line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(key, argname, starred, ','.join(prots), str(len(prots)),','.join(set(orglist)),num_isolates)
		outfile.write(line)
