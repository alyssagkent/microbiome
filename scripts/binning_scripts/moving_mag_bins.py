from shutil import copyfile
#with open file
#B314-1.004_sub.contigs k__Bacteria 0.85 0.00 0.00 240526 75.0 562 species 7 18.75 k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia_coli;
#if it is
#wrk

wrk = '/workdir/users/agk85/CDC2/bins'
bin_handle = '/workdir/users/agk85/CDC2/bins/all_das_checkm_kraken.txt'
with open(bin_handle) as infile:
	for line in infile:
		fields = line.strip().split('\t')
		name = fields[0]
		completeness = float(fields[2])
		contamination = float(fields[3])
		genome_size = float(fields[5])
		kraken_completeness = float(fields[6])
		kraken_contamination = float(fields[10])
		kraken_taxonomy = fields[11]
		kraken_rank = fields[8]
		#sort the file into different things
		#start with the checkm stuff first
		src = wrk + '/all_bins/' + name + '.fa'
		#hq bins --- checkm >90% complete, <5% contamination
		if (completeness >90 and contamination <5):
			#cp file to hq folder
			dst = wrk + '/hq_bins/' + name + '.fa'
			copyfile(src, dst)			
		#mq bins --- checkm >50% complete, <10% contamination
		if (completeness >50 and contamination <10):
			#cp file to mq folder
			dst = wrk + '/mq_bins/' + name  + '.fa'
			copyfile(src, dst)			
		#all quality bins --- <10% contamination
		if (contamination <10):
			#cp file to quality folder
			dst = wrk + '/quality_bins/' + name  + '.fa'
			copyfile(src, dst)	
		#hq species ---- checkm > 90% complete, <5% contamination, kraken_rank == species
		if (completeness >90 and contamination <5 and kraken_rank == 'species'):
			species = kraken_taxonomy.split('s__')[1].split(';')[0]
			species = species.replace('/', '?')
			dst = wrk + '/hq_species_bins/' + species + '_' + name + '.fa'
			copyfile(src, dst)  
		
		
