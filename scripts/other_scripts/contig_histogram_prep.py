
# WRK=/workdir/users/agk85/CDC2
# NAME=B314-2
# python contig_histogram_prep.py \
# -o $WRK/bins/${NAME}_contig_taxa_das_2_fgs \
# -i $WRK/bins/${NAME}_das_2_contigs_v_bins.txt \
# -c $WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
#files should come out as arg_orgcounts_${NAME}_species1.txt or species0.txt if not including hic

import argparse
from argparse import RawDescriptionHelpFormatter

#def getOptions():
description="""Programs reads through contig-bin file and aggregates taxa at particular levels"""

parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i','--infile', dest="inhandle", action='store', required=True,  help='Infile [REQUIRED]', metavar="INFILE")
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED]', metavar="OUTFILE")
parser.add_argument('-c','--combo', dest="combohandle", action='store', required=True, help='Combo table [REQUIRED]', metavar="COMBOTABLE")

args = parser.parse_args()


percent_thresh = 0.1
contact_thresh = 2

levels = ['f__','g__', 's__']

contig_taxa_count = {}
#big dict
for level in levels:
	contig_taxa_count[level]={}


contigs = []
with open(args.combohandle) as combofile:
	header = combofile.readline()
	for line in combofile:
		contig = line.split('\t')[0]
		contigs.append(contig)


threshes = ['2+contacts','10+percent']
bintypes = ['quality','anybin','bin+contig']
#initiation
contig_taxa = {}
for thresh in threshes:
	contig_taxa[thresh] = {}
	for bintype in bintypes:
		contig_taxa[thresh][bintype]={}
		for contig in contigs:
			contig_taxa[thresh][bintype][contig] = []
			

#iterate through file gathering bin taxa
#count   sample  contig  bin     bin_length      quality association_type        total_count     trans_count     norm_l  norm_rf contig_taxonomy kraken_taxonomy gtdb_taxonomy   arg_prsence     arg_clusters    arg_genes       mge_presence    mge_clusters    mge_genes
#cluster_resident, cluster_hic, contig_resident, contig_hic
#NA, BAD, LQ, MQ, HQ

total_arg_contigs = 0
with open(args.inhandle) as infile:
	header = infile.readline()
	wonkkyline = infile.readline()
	for line in infile:
		arg = line.split('\t')[14]
		if arg == '1':
			total_arg_contigs += 1
			qualitybin = 0
			anybin = 1
			bin_or_contig = 0
			contacts_2 = 1
			percent_10 = 0
			quality = line.split('\t')[5]
			try:
				prop = float(line.split('\t')[8])/float(line.split('\t')[7])
			except ZeroDivisionError:
				prop = 0
			taxonomy = line.split('\t')[12]
			contig = line.split('\t')[2]
			contacts = float(line.split('\t')[8])
			if prop >= percent_thresh:
				if (quality == 'HQ' or quality == 'MQ'):
					contig_taxa['10+percent']['quality'][contig].append(taxonomy)
				if (quality == 'BAD' or quality == 'HQ' or quality == 'MQ' or quality == 'LQ'):
					contig_taxa['10+percent']['anybin'][contig].append(taxonomy)
				contig_taxa['10+percent']['bin+contig'][contig].append(taxonomy)				
			if contacts >= contact_thresh:
				if (quality == 'HQ' or quality == 'MQ'):
					contig_taxa['2+contacts']['quality'][contig].append(taxonomy)
				if (quality == 'BAD' or quality == 'HQ' or quality == 'MQ' or quality == 'LQ'):
					contig_taxa['2+contacts']['anybin'][contig].append(taxonomy)
				contig_taxa['2+contacts']['bin+contig'][contig].append(taxonomy)


sample = args.combohandle.split('_master_scf_table.txt')[0].split('/')[-1]
header = 'Count\tSample\tThresh\tBintype\tLevel\tContig\tTaxaCount\tTotal_arg_contigs\tTaxaset\n'
count = 0
with open(args.outhandle,'w') as outfile:
	outfile.write(header)
	for thresh in threshes:
		for bintype in bintypes:
			for level in levels:
				for contig in contigs:
						taxalist = []
						for taxonomy in contig_taxa[thresh][bintype][contig]:
							if taxonomy != 'NA' and taxonomy != '.':
								taxa = taxonomy.split(level)[1].split(';')[0]
								uptotaxa = taxonomy.split(level)[0] + level + taxonomy.split(level)[1].split(';')[0]
								if taxa != '':
									taxalist.append(uptotaxa)
						taxaset = set(taxalist)
						taxasetlen = len(taxaset)
						if taxasetlen > 0:
							outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(str(count), sample, thresh, bintype, level, contig, str(taxasetlen), str(total_arg_contigs), ','.join(taxalist)))
						count += 1