
#get_contig_arg_plasmid_phage_tpn_counts.py
#Cluster Pfam_annotation MGE_type        MGE_combo_type  Plasmid Phage   Transposon      ICE     IME     Other   ARG     Cluster2        Rep_gene        Genes   Gene_count      Pfam_ids


#contig dictionary
contig_dict={}

def update_contig_counts(contig_dict, contig, gene_info):
	try:
		contig_dict[contig]= tuple(x + y for x,y in zip(contig_dict[contig], gene_info))
	except KeyError:
		contig_dict[contig] = (0,0,0,0,0,0)

#open the contig file
with open('~/agk/CDC2/mobile/metagenomes/old_mge_99_nr.fna.clstr.tbl.annot') as infile:
	header = infile.readline()
	for line in infile:
		fields= line.split('\t')
		genes = line.split('\t')[13].split(',')
		gene_info = (int(fields[4]),int(fields[5]),int(fields[6]),int(fields[7]),int(fields[8]),int(fields[10]))
		for gene in genes:
			contig = '_'.join(gene.split('_')[0:-1])
			update_contig_counts(contig_dict,contig, gene_info)


with open('~/agk/CDC2/mobile/metagenomes/contig_mge_counts.txt','w') as outfile:
	outfile.write('Contig\tPlasmid\tPhage\tTransposon\tICE\tIME\tARG\n')
	for contig in contig_dict.keys():
		gene_info = contig_dict[contig]
		print(contig)
		print(gene_info)
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(contig, str(gene_info[0]),str(gene_info[1]), str(gene_info[2]),str(gene_info[3]),str(gene_info[4]),str(gene_info[5])))
