#make a more parse-able table from cd-hit output
#real specific right now for the mge's
import sys
import glob

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

pfamiddesc = {}
pfaminfo = open('/workdir/users/agk85/CDC/tables/metagenomes3/pfam_acc_desc.txt')
for line in pfaminfo:
	pfamid = line.split('\t')[0].split('ACC   ')[1]
	desc = line.split('\t')[1].split('DESC  ')[1].strip()
	pfamiddesc[pfamid] = desc

pfamdict = {}
pfamiddict = {}
pfamdescdict = {}
pfampaths = glob.glob('/workdir/users/agk85/CDC/annotation/metagenomes3/*_pfam.txt')
for pfamfile in pfampaths:
	with open(pfamfile) as p:
		for line in p:
			pfamid = line.split('\t')[3]
			desc = pfamiddesc[pfamid]
			pfamdescdict[line.split('\t')[0]] = desc


#make a tab delimited file in this format
name = 'cluster\targname\tstarredid\tids\tnumber_of_proteins\torgs\tdistinct_isolate_species\n'
with open(outhandle, 'w') as outfile:
	for key in clusterprotdict.keys():
		prots = clusterprotdict[key]
		starred = clusternum_map[key]
		pfamname = ''
		orgs = ''
		orglist = []
		mgmlist = []
		for prot in clusterprotdict[key]:
			if 'NODE' not in prot:
				mgmlist.append(prot.split('_')[0])
				try:
					pfamname = pfamname + ',' + pfamdescdict[starred].replace(' ', '_')
				except:
					pfamname = pfamname
			if 'NODE' in prot:
				isolate = prot.split('_')[0]
				org = isolate_orgs[isolate]
				orgs = orgs + ',' + org
				orglist.append(org)
		orgclean = ','.join(set(orglist)).replace(' ', '_')
		if orgclean == '':
			orgclean = 'NA'
		if pfamname == '':
			pfamname = 'NA'
		num_isolates = str(len(set(orglist)))	
		num_mgms = str(len(set(mgmlist)))
		line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(key, pfamname, starred, ','.join(prots), str(len(prots)),orgclean,num_isolates,num_mgms)
		outfile.write(line)
