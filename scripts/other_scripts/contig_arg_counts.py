
#contig_arg_counts.py
#4973    B331-3_scaffold_414_28  B331-3_scaffold_414_28,B331-4_scaffold_458_18   2


#contig dictionary
contig_dict={}

def update_contig_counts(contig_dict, contig, gene_info):
	try:
		contig_dict[contig]= tuple(x + y for x,y in zip(contig_dict[contig], gene_info))
	except KeyError:
		contig_dict[contig] = (0)

#open the contig file
with open('~/agk/CDC2/arg_v_org/metagenomes/old/args_99_nr.fna.clstr.tbl') as infile:
	header = infile.readline()
	for line in infile:
		fields= line.split('\t')
		genes = line.split('\t')[2].split(',')
		#gene_info = (1)
		for gene in genes:
			contig = '_'.join(gene.split('_')[0:-1])
			try:
				contig_dict[contig] +=1
			except KeyError:
				contig_dict[contig] = 1
			#update_contig_counts(contig_dict,contig, gene_info)


with open('~/agk/CDC2/arg_v_org/metagenomes/old/contig_arg_counts.txt','w') as outfile:
	outfile.write('Sample\tContig\tARG\n')
	for contig in contig_dict.keys():
		gene_info = contig_dict[contig]
		print(contig)
		print(gene_info)
		sample = contig.split('_')[0]
		outfile.write('{0}\t{1}\t{2}\n'.format(sample, contig, str(gene_info)))
