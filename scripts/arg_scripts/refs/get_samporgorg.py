

def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())



#dictionary
#first initiate the dictionary with the samples and the clusters

#first the samples
samples = []
sampargorg = {}
samporgorg = {}
with open('/workdir/users/agk85/CDC/HicDesign.txt') as infile:
	for line in infile:
		sample = line.strip()
		sampargorg[sample] = {}
		samporgorg[sample] = {}
		samples.append(sample)

#get all the relevant clusters
clusters = []
with open("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_org_hic_cliques_95_98_1_2.tbl") as cliques:
	header = cliques.readline().strip().split('\t')
	#Count   Cluster ARG_name        Top_ARG Sample  k__Bacteria; p__; c__; o__; f__; g__; s__;
	orgs = header[5:]
	for line in cliques:
		cluster = line.split('\t')[1]
		clusters.append(cluster)

#make them a set and sort them
clusterset = list(set(clusters))
clusterset.sort(natcasecmp)

for sample in samples:
	for cluster in clusterset:
		sampargorg[sample][cluster] = []


count = 0
#now go through and add in all of the organisms
with open("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_org_hic_cliques_95_98_1_2.tbl") as cliques:
	header = cliques.readline().strip().split('\t')
	#Count   Cluster ARG_name        Top_ARG Sample  k__Bacteria; p__; c__; o__; f__; g__; s__;
	orgs = header[5:]
	for line in cliques:
		orghits = line.strip().split('\t')[5:]
		for i in range(len(orghits)):
			orghit = orghits[i]
			if orghit == '1':
				count = count + 1
				org = orgs[i]
				sample  = line.split('\t')[4]
				cluster = line.split('\t')[1]
				sampargorg[sample][cluster].append(org)



tried = []
#now double loop through organisms and then loop through sample and org and if it has both organisms then add one to it
for org1 in orgs:
	for org2 in orgs:
		if org2+org1 not in tried:
			tried.append(org1+org2)
			for sample in samples:
				samporgorg[sample][org1+org2] = 0
				for cluster in clusterset:
					if org1 in sampargorg[sample][cluster]:
						if org2 in sampargorg[sample][cluster]:
							samporgorg[sample][org1+org2] = samporgorg[sample][org1+org2] + 1

tried = []
header = ''
for org1 in orgs:
	for org2 in orgs:
		if org1 != org2:
			if org2+org1 not in tried:
				tried.append(org1+org2)
				header = header + org1 + "+" + org2 + '\t' 

header = header + '\n'

tried = []
with open('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/sample_vs_orgargorg_cliques_95_98_1_2.tbl','w') as outfile:
	outfile.write(header)
	for sample in samples:
		print(sample)
		outfile.write(sample)
		for org1 in orgs:
			for org2 in orgs:
				if org1 != org2:
					if org2+org1 not in tried:
						tried.append(org1+org2)
						outfile.write('\t{0}'.format(str(samporgorg[sample][org1+org2])))
		outfile.write('\n')
	
						
