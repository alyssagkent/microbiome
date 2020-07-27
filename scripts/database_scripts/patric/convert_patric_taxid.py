#converting file to have the taxid in the beginning

#/home/agk85/wd/refdbs/PATRIC
#not sure how i got to the genome_ids.txt---presumably it comes from the PATRIC database

import sys
gene_to_genomeid = {}
with open('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/patric_comparisons/genome_ids.txt') as infile:
	for line in infile:
		gene = line.split(' ')[0].split('>')[1]
		genome = line.split('.fna')[0]
		gene_to_genomeid[gene] = genome
		

genome_to_taxid= {}
with open('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/patric_comparisons/genome_metadata') as infile:
	for line in infile:
		genomeid = line.split('\t')[0]
		taxid = line.split('\t')[3]
		genome_to_taxid[genomeid] = taxid


taxid_to_lineage= {}
with open('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/patric_comparisons/genome_lineage') as infile:
	header = infile.readline()
	for line in infile:
		lineage = line.strip().split('\t')[11].split(';')
		taxid = line.split('\t')[2]
		try:
			newlineage = 'k__{0}; group__{1}; p__{2}; c__{3}; o__{4}; f__{5}; g__{6}; s__{7};'.format(lineage[1],lineage[2],lineage[3],lineage[4],lineage[5],lineage[6],lineage[7],lineage[8])
			newlineage = lineage
			taxid_to_lineage[taxid] = newlineage
		except IndexError:
			taxid_to_lineage[taxid] = lineage


#taxids = ['702953', '2211115', '702969', '702959', '2250596']
#for taxid in taxids:
#	print taxid, taxid_to_lineage[taxid]	
	
	

#ARGs_Hi-C_BLASTn_all_1e-100.txt
inhandle = sys.argv[1]
outhandle = sys.argv[2]
with open(inhandle) as infile:
	with open(outhandle,'w') as outfile:
		header = infile.readline()
		for line in infile:
			ref = line.split('\t')[1]
			genomeid = gene_to_genomeid[ref]
			taxid = genome_to_taxid[genomeid]
			outfile.write(taxid + '\t' + line)

