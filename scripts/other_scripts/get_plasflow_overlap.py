#USAGE python get_plasflow_overlap.py /workdir/users/agk85/press2/combo_tables/metagenomes/ProxiMeta-1_master_scf_table.txt /workdir/users/agk85/press2/combo_tables/metagenomes/ProxiMeta-1_plasflow_overlap.txt /workdir/users/agk85/press2/combo_tables/metagenomes/ProxiMeta-1_plasmid_lengths.txt ProxiMeta-1

import sys
inhandle = sys.argv[1]
outhandle = sys.argv[2]
outhandle2 = sys.argv[3]
sample = sys.argv[4]
keys = ['a','b','c','d','e','f','g','h']
total_dict ={'a':0, 'b':0, 'c':0, 'd':0, 'e':0, 'f':0, 'g':0, 'h':0}
shared_dict={'a':0, 'b':0, 'c':0, 'd':0, 'e':0, 'f':0, 'g':0, 'h':0}
total_len_dict={'a':0, 'b':0, 'c':0, 'd':0, 'e':0, 'f':0, 'g':0, 'h':0}

with open(inhandle) as infile:
	header = infile.readline()
	for line in infile:
		contig_length = float(line.split('length_')[1].split('\t')[0])
		if '1plasmid' in line:
			total_dict['a'] +=1
			total_len_dict['a'] +=contig_length
			if '8plasmid' in line:
				shared_dict['a'] +=1
		if '2plasmid' in line:
			total_dict['b'] +=1
			total_len_dict['b'] +=contig_length
			if '8plasmid' in line:
				shared_dict['b'] +=1
		if '3plasmid' in line:
			total_dict['c'] +=1
			total_len_dict['c'] +=contig_length
			if '8plasmid' in line:
				shared_dict['c'] +=1
		if '4plasmid' in line:
			total_dict['d'] +=1
			total_len_dict['d'] +=contig_length
			if '8plasmid' in line:
				shared_dict['d'] +=1
		if '5plasmid' in line:
			total_dict['e'] +=1
			total_len_dict['e'] +=contig_length
			if '8plasmid' in line:
				shared_dict['e'] +=1
		if '6plasmid' in line:
			total_dict['f'] +=1
			total_len_dict['f'] +=contig_length
			if '8plasmid' in line:
				shared_dict['f'] +=1
		if '7plasmid' in line:
			total_dict['g'] +=1
			total_len_dict['g'] +=contig_length
			if '8plasmid' in line:
				shared_dict['g'] +=1
		if '8plasmid' in line:
			total_dict['h'] +=1
			total_len_dict['h'] +=contig_length
			if '8plasmid' in line:
				shared_dict['h'] +=1


header = 'Sample\tPlasmid_Finder\tFull_Plasmids\tRelaxase\tPlasmid_pfam\tImme_plasmid\tAclame_plasmid\tPfams_plasmid\tPlasflow\n'
with open(outhandle,'w') as outfile:
	outfile.write(header)
	outfile.write(sample)
	for key in keys:
		outfile.write('\t=' + str(shared_dict[key]) + '/' + str(total_dict[key]))
	outfile.write('\n')

with open(outhandle2,'w') as outfile:
	outfile.write(header)
	outfile.write(sample)
	for key in keys:
		outfile.write('\t=' + str(total_len_dict[key]) + '/' + str(total_dict[key]))
	outfile.write('\n')	
