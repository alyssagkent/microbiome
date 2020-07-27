#this program appends the number of contigs within each bin that gets hit 1 otherwise

#bindict
bindict = {}
with open('/workdir/users/agk85/CDC2/das/all_DASTool_scaffolds2bin.txt') as binfile:
	for line in binfile:
		contig,binid = line.strip().split('\t')
		try:
			bindict[binid].append(contig)
		except KeyError:
			bindict[binid] = [contig]

partners = {}
with open('/workdir/users/agk85/CDC2/hicpro/CDC+healthy_trans_primary_ncol_0_noeuks.txt') as transfile:
	for line in transfile:
		contig1,contig2,count=line.strip().split('\t')
		try:
			partners[contig1].append(contig2)
		except KeyError:
			partners[contig1]=[contig2]
		try:
			partners[contig2].append(contig1)
		except KeyError:
			partners[contig2]=[contig1]


def number_contigs_hit(contig, binid):
	#access some dicts
	num_contigs = 0
	binidmembers = bindict[binid]
	try:
		partnerlist = list(set(partners[contig]))
		for partner in partnerlist:
			if partner in binidmembers:
				num_contigs+=1 	
	except KeyError:
		a = 0
	return num_contigs


with open('/workdir/users/agk85/CDC2/bins/all_das_1_contigs_v_bins_all.txt') as infile, open('/workdir/users/agk85/CDC2/bins/all_das_1_contigs_v_bins_contignum.txt','w') as outfile:
	outfile.write(infile.readline().strip()+ '\tnumber_contigs_hit\n')
	for line in infile:
		contig = line.split('\t')[2]
		binid = line.split('\t')[3]
		association_type = line.split('\t')[6]
		if association_type == 'cluster_hic' or association_type == 'cluster_resident':
			bincontigshit = number_contigs_hit(contig,binid)
			outfile.write(line.strip() + '\t' + str(bincontigshit) + '\n')
		else:
			if contig == binid:
				outfile.write(line.strip() + '\t0\n')
			else:
				outfile.write(line.strip() + '\t1\n')



