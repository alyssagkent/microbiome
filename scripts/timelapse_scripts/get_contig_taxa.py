#lame script to get whether or not the taxonomy I think the bin has is represented by the 

contigtaxadict = {}
with open('/workdir/users/agk85/CDC2/kraken/all.kraken.taxonomy.txt') as infile:
	for line in infile:
		contig,ktaxonomy = line.strip().split('\t')
		contigtaxadict[contig] = ktaxonomy


markerlist = ['d','p','c','o','f','g','s']
base_marker = {'d':'','p':'','c':'','o':'','f':'','g':'','s':''}
def convert_taxonomy(ktaxonomy):
	taxonomy_fields = ktaxonomy.split('|')
	for field in taxonomy_fields:
			if field != 'root':
					marker = field.split('__')[0]
					taxa = field.split('__')[1]
					if marker in markerlist:
							base_marker[marker] = taxa
	parsedtaxonomy = ''
	for marker in markerlist:
			#standard marker is 'k' for kingdom instead of domain, so switching that
			if marker == 'd':
					parsedtaxonomy = parsedtaxonomy + 'k__' + base_marker[marker]
			else:
					parsedtaxonomy = parsedtaxonomy + marker + '__' + base_marker[marker]
			#my taxonomy format stupidly has semicolon then spaces except after the species
			if marker != 's':
					parsedtaxonomy = parsedtaxonomy + '; '
			else:
					parsedtaxonomy = parsedtaxonomy + ';'
	return parsedtaxonomy


trans = {}
with open('/workdir/users/agk85/CDC2/hicpro/CDC+healthy_trans_primary_ncol_0_noeuks.txt') as infile:
	for line in infile:
		contig1,contig2,count= line.strip().split('\t')
		try:
			trans[contig1].append(contig2)
		except KeyError:
			trans[contig1] = [contig2]
		try:
			trans[contig2].append(contig1)
		except KeyError:
			trans[contig2] = [contig1]

def get_presence_rel_irr(contiglist,gtaxonomy):
	taxapresent = 0
	relevant_combos=[]
	irrelevant_combos = []
	for genecontig in contiglist:
		try:
			contigs = trans[genecontig]
			for taxacontig in contigs:
				ktaxonomy = contigtaxadict[taxacontig]
				ctaxonomy = convert_taxonomy(ktaxonomy)
				if gtaxonomy == ctaxonomy:
					taxapresent = 1
					relevant_combos.append(taxacontig + ':' + ctaxonomy)
				irrelevant_combos.append(taxacontig + ':' + ctaxonomy)
		except KeyError:
			a = 1
			#there are no hic reads
			relevant_combos.append('no_trans')
			irrelevant_combos.append('no_trans')
	rel_combos = ','.join(relevant_combos)
	irr_combos = ','.join(irrelevant_combos)
	if rel_combos == '':
		rel_combos = 'no_taxa'
	if irr_combos == '':
		irr_combos = 'no_taxa'
	return str(taxapresent),rel_combos, irr_combos

inhandle = "/workdir/users/agk85/CDC2/bins/timelapse/timelapse_mge_org_2_alltaxalevels_gained_filteredindex.txt"
outhandle= "/workdir/users/agk85/CDC2/bins/timelapse/timelapse_mge_org_2_alltaxalevels_gained_filteredindex_krakencontig.txt"
with open(inhandle) as infile,open(outhandle,'w') as outfile:
	header = infile.readline()
	header = header + '\tt2_contigtaxa_presence\tt2_contigtaxa_rel\tt2_contigtaxa_irr\n'
	outfile.write(header)
	for line in infile:
		count,connection_count,patient1,t1,t2,level,taxa1,taxa1_fulltaxa,ho,taxa1_t1_abund,taxa1_t2_abund,taxa2,fulltaxa2,taxa2_t1_abund,taxa2_t2_abund, gene,genename,genecontigs_t2,t1rpkm, t2rpkm,trans_taxa2_t1,trans_taxa2_t2,num_nont1,other_tps,contig_taxa_connection,contigfulltaxa1_t1,t1ll,trans_taxa1_t1,trans_taxa1_t2,num_bincontigs_t1, num_bincontigs_t2,binids_t1,binids_t2,filterstatus,filter_reason,higher_order,higher_level_t1,higher_level_t2,trans_taxa2_t1_morethan_0,trans_taxa2_t2_lessthan_2,t1contigpres,trans__taxa1_t1_equals_0,genetype= line.split('\t')
		taxapresent = 0
		relevant_combos=[]
		irrelevant_combos = []
		#presence_t1,rel_t1,irr_t1=get_presence_rel_irr(genecontigs_t1)
		presence_t2,rel_t2,irr_t2=get_presence_rel_irr([genecontigs_t2],fulltaxa2)
		#get the contigs that map to genecontig #i don't think you need to check beforehand
		#for each, ask if the kraken contig taxonomy agrees with the bin taxonomy ---some obviously won't but this 
		#if it does, then add it to the list for this contig
		outfile.write('{0}\t{1}\t{2}\t{3}\n'.format(line.strip(), presence_t2,rel_t2,irr_t2))

