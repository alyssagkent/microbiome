#!/usr/local/bin/python
description="""This script is to merge our bin taxonomies with metaphlan abundances....as best as we can do"""
# print(description)
#######################################

#load bintables
#load metaphlan
#merge on every taxa level taking care to have the right formatting
# if you can't find an exact one, then ummmmmmmmmmmmmm write NA 

# get a way to fill in the rest of the taxonomy depending on what level it is
therest = {'k':' p__; c__; o__; f__; g__; s__;','p':' c__; o__; f__; g__; s__;','c':' o__; f__; g__; s__;','o':' f__; g__; s__;','f':' g__; s__;','g':' s__;','s':''}

metaphlandict = {}
with open('/workdir/users/agk85/CDC2/metaphlan/cdc/mgm/CDC_mgm_metaphlan.txt') as mfile:
	header = mfile.readline()
	samples = header.strip().split('\t')
	sid = samples.pop(0)
	for sample in samples:
		metaphlandict[sample] = {}
	
	for line in mfile:
		abunds = line.strip().split('\t')
		taxonomy= abunds.pop(0)
		mytaxonomy = taxonomy.replace('|','; ') + ';'
		mytaxonomy = mytaxonomy.replace('Enterobacteriales','Enterobacterales')
		mytaxonomy = mytaxonomy.replace('Eubacterium_hallii','Anaerobutyricum_hallii')
		mytaxonomy = mytaxonomy.replace('c__Actinobacteria','c__Coriobacteriia')
		mytaxonomy = mytaxonomy.replace('Coriobacteriales','Eggerthellales')
		mytaxonomy = mytaxonomy.replace('Coriobacteriaceae','Eggerthellaceae')
		mytaxonomy = mytaxonomy.replace('o__Selenomonadales; f__Veillonellaceae','o__Veillonellales; f__Veillonellaceae')
		level = mytaxonomy.split('; ')[-1].split('__')[0]
		shorttaxonomy = level + '__' + mytaxonomy.split(level + '__')[1].split(';')[0]
		if level != 't':
			fulltaxa = mytaxonomy + therest[level]
			for sample, abund in zip(samples, abunds):
				metaphlandict[sample][shorttaxonomy] = abund

bindict = {}
for sample in samples:
	bindict[sample] = {}
levels = ['k__','p__','c__','o__','f__','g__','s__']
levelnums = {'k__':1,'p__':2,'c__':3,'o__':4,'f__':5,'g__':6,'s__':7}
binhandle = '/workdir/users/agk85/CDC2/das/all_bintables.txt'
outhandle = '/workdir/users/agk85/CDC2/das/all_bintables_metaphlan.txt'
header = 'binid\tpatient\tsample\tquality\ttaxonomy\tk__\tp__\tc__\to__\tf__\tg__\ts__\n'
with open(binhandle) as binfile,open(outhandle,'w') as outfile:
	for line in binfile:
		binid,patient,sample,quality,taxonomy = line.strip().split('\t')
		outfile.write(line.strip())
		final_valid_binlevel = 0
		final_valid_metaphlanlevel = 0
		for level in levels:
			try:
				taxa=taxonomy.split(level)[1].split(';')[0]
				if (taxa != '' and taxa != '.'):
					final_valid_binlevel = levelnums[level]
				taxon = level + taxa
				abundance = metaphlandict[sample][taxon]
				final_valid_metaphlanlevel = levelnums[level]
			except :
				abundance = 'NA'
			outfile.write('\t' + abundance)
		outfile.write('\t' + str(final_valid_binlevel) + '\t' + str(final_valid_metaphlanlevel) + '\t' + str(final_valid_binlevel-final_valid_metaphlanlevel))
		outfile.write('\n')


			


