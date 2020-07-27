#bin_cleanliness.py

#WRK=/workdir/users/agk85/CDC2
#NAME=B357-3
# python bin_cleanliness.py \
# -b $WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt \
# -o $WRK/bins/${NAME}_bin_cleanliness.txt \
# -l $WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_allValidPairs \
# -t $WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt \
# -k $WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt \
# -c $WRK/das/${NAME}/checkm_lineage/${NAME}.stats

import argparse
from argparse import RawDescriptionHelpFormatter
from copy import deepcopy

#def getOptions():
description="""
#This program will go through all of the reads for a sample and break into compartments with and without mges
"""
parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="INFILE")
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED]', metavar="OUTFILE")
parser.add_argument('-l','--hic', dest="hichandle", action='store', required=True, help='HI-C file [REQUIRED]', metavar="HIC_FILE")
parser.add_argument('-t','--table', dest="combotable", action='store', required=True, help='ComboTable file [REQUIRED]', metavar="COMBO_FILE")
parser.add_argument('-c','--checkm', dest="checkm", action='store', required=True, help='checkm_info [REQUIRED]', metavar="CHECKM")
#parser.add_argument('-m','--min', dest="min", action='store', type=int, required=True, help='Min contacts [REQUIRED]', metavar="MIN_PROP")
parser.add_argument('-k','--kraken', dest="kraken_bins", action='store', required=False, help='Kraken cluster file', metavar="KRAKEN")

args = parser.parse_args()
sample = args.kraken_bins.split('/')[-1].split('_all_kraken_weighted.txt')[0]

contigs = {}
mgedict = {}
noncluster_contigs  =  []
#get all the contigs and if its a mobile contig
with open(args.combotable) as infile:
	header = infile.readline()
	for line in infile:
		contig = line.split('\t')[0]
		noncluster_contigs.append(contig)
		mge = line.split('\t')[-2]
		contigs[contig] = 1
		mgedict[contig] = mge

#checkm quality
qualitydict = {}
cluster_length = {}
with open(args.checkm) as checkm:
	for line in checkm:
		bin = line.split('\t')[0].split('.contigs')[0]
		completion = float(line.split('\t')[2])
		contamination = float(line.split('\t')[3])
		length = line.strip().split('\t')[5]
		#print(completion, contamination)
		quality = 'BAD'
		if (completion > 90 and contamination < 5):
			quality = 'HQ'
		elif (completion >= 50 and contamination < 10):
			quality = 'MQ'
		elif (completion < 50 and contamination < 10):
			quality = 'LQ'	
		qualitydict[bin] = quality
		cluster_length[bin] = length


#start with full list and remove iteratively as you find them in the bins
#do not include bad bins
####BIN STUFF######
#create a dictionary of the contigs and bins
bincontigsdict = {} #put in bin and get contigs
bindict = {} #put in contig and get bin (meaning that some will not have a bin)
with open(args.binhandle) as binfile:
	for line in binfile:
		bin = line.strip().split('\t')[1]
		contig = line.split('\t')[0]
		if qualitydict[bin] != 'BAD':
			bindict[contig] = bin
			#get rid of this contig (because it has a bin)
			try:
				noncluster_contigs.remove(contig)
			except:
				a =1 #it's a kraken euk!

#print(bindict)
for contig in noncluster_contigs:
	bindict[contig] = 'binless'


#kraken bin taxonomies	
kraken_bin_taxonomy = {}
with open(args.kraken_bins) as krakenfile:
	for line in krakenfile:
		bin = line.split('.contigs.fa.report.txt.besttaxid')[0]
		taxonomy = line.strip().split('\t')[7]
		kraken_bin_taxonomy[bin] = taxonomy
		
#get hic info
#A00257:145:HF5YNDRXX:1:2168:15745:24439 B314-1_1|phage|127443|142740_20775      101     -       B314-1_scaffold_1       127362  +       16000   HIC_B314-1_1|phage|127443|142740_20775_1        HIC_B314-1_scaffold_1_306       42      42
#ok so there are some intricacies
#imagine read1 and read2
#they could hit
#a) same contig, same bin  -- this increases intra_bin_cis
#b) same contig, binless -- this increases intra_binless_cis

#e) different contigs, same bin -- this increases intra_bin_trans-
#f) different contigs, both binless -- this increases intra_binless_trans-
#g) different contigs, different bins -- this increases inter_bin_trans
#h) different contigs, one bin one binless -- this increases inter_binless_trans

intra_bin_cis = 0
intra_binless_cis =0
intra_bin_trans = 0
intra_binless_trans = 0
inter_bin_trans = 0
inter_binless_trans = 0
header = 'sample\tmobility\ttotal_reads\ttotal_trans\tintra_bin_cis\tintra_binless_cis\tintra_bin_trans\tintra_binless_trans\tinter_bin_trans\tinter_binless_trans\tintra_bin_cis_prop\tintra_binless_cis_prop\tintra_bin_trans_prop\tintra_binless_trans_prop\tinter_bin_trans_prop\tinter_binless_trans_prop\n'
#create bin table with hic data (long format)
with open(args.outhandle, 'w') as outfile:
	outfile.write(header)
	with open(args.hichandle) as hicfile:
		for line in hicfile:
			#Get the contigs linking
			contig1 = line.split('\t')[1]
			contig2 = line.split('\t')[4]
			try:
				a= contigs[contig1]
				b = contigs[contig2]
				#deal with phage
				#B314-1_1|phage|127443|142740_20775 -> B314-1_scaffold_1
				if 'phage' in contig1:
					parent = contig1.split('_')[0] + '_scaffold_' + contig1.split('_')[1].split('|')[0]
					if parent == contig2:
						contig1 = contig2
				if 'phage' in contig2:
					parent = contig2.split('_')[0] + '_scaffold_' + contig2.split('_')[1].split('|')[0]
					if parent == contig1:
						contig2 = contig1		
				#split on if it is same or different contig
				if contig1 == contig2: #they are cis
					bin1 = bindict[contig1]
					if bin1 == 'binless':
						intra_binless_cis += 1
					else:
						intra_bin_cis += 1
				else: #they are trans
					bin1 = bindict[contig1]
					bin2 = bindict[contig2]
					if bin1 == bin2:
						if bin1 == 'binless':
							intra_binless_trans += 1
						else:
							intra_bin_trans += 1
					else: #different bins
						if ((bin1=='binless') or (bin2 == 'binless')):
							inter_binless_trans+=1
						else:
							inter_bin_trans+=1
			except KeyError:
				print('one is a euk', contig1, contig2)
		mobility = 'withmge'
		total_reads = intra_bin_cis+intra_binless_cis+intra_bin_trans+intra_binless_trans+inter_bin_trans+inter_binless_trans
		total_trans = intra_bin_trans+intra_binless_trans+inter_bin_trans+inter_binless_trans
                intra_bin_cis_prop = intra_bin_cis/float(total_reads)
                intra_binless_cis_prop = intra_binless_cis_prop/float(total_reads)
                intra_bin_trans_prop = intra_bin_trans/float(total_reads)
                intra_binless_trans_prop= intra_binless_trans/float(total_reads)
                inter_bin_trans_prop= inter_bin_trans/float(total_reads)
                inter_binless_trans_prop= inter_binless_trans/float(total_reads)
                outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n'.format(sample, mobility, total_reads, total_trans, intra_bin_cis, intra_binless_cis, intra_bin_trans, intra_binless_trans, inter_bin_trans, inter_binless_trans, intra_bin_cis_prop, intra_binless_cis_prop, intra_bin_trans_prop, intra_binless_trans_prop, inter_bin_trans_prop, inter_binless_trans_prop))
	#do it again but without mges
	with open(args.hichandle) as hicfile:
		intra_bin_cis = 0
		intra_binless_cis =0
		intra_bin_trans = 0
		intra_binless_trans = 0
		inter_bin_trans = 0
		inter_binless_trans = 0
		for line in hicfile:
			#Get the contigs linking
			contig1 = line.split('\t')[1]
			contig2 = line.split('\t')[4]
			try:
				a = contigs[contig1]
				b = contigs[contig2]				
				mge1 = mgedict[contig1]
				mge2 = mgedict[contig2]
				if ((mge1!='mge') and (mge2!='mge')):
					#deal with phage
					#B314-1_1|phage|127443|142740_20775 -> B314-1_scaffold_1
					if 'phage' in contig1:
						parent = contig1.split('_')[0] + '_scaffold_' + contig1.split('_')[1].split('|')[0]
						if parent == contig2:
							contig1 = contig2
					if 'phage' in contig2:
						parent = contig2.split('_')[0] + '_scaffold_' + contig2.split('_')[1].split('|')[0]
						if parent == contig1:
							contig2 = contig1		
					#split on if it is same or different contig
					if contig1 == contig2: #they are cis
						bin1 = bindict[contig1]
						if bin1 == 'binless':
							intra_binless_cis += 1
						else:
							intra_bin_cis += 1
					else: #they are trans
						bin1 = bindict[contig1]
						bin2 = bindict[contig2]
						if bin1 == bin2:
							if bin1 == 'binless':
								intra_binless_trans += 1
							else:
								intra_bin_trans += 1
						else: #different bins
							if ((bin1=='binless') or (bin2 == 'binless')):
								inter_binless_trans+=1
							else:
								inter_bin_trans+=1
			except KeyError:
				print('one is a euk', contig1, contig2)
		mobility = 'nomge'
		total_reads = intra_bin_cis+intra_binless_cis+intra_bin_trans+intra_binless_trans+inter_bin_trans+inter_binless_trans
		total_trans = intra_bin_trans+intra_binless_trans+inter_bin_trans+inter_binless_trans
		intra_bin_cis_prop = intra_bin_cis/float(total_reads)
		intra_binless_cis_prop = intra_binless_cis_prop/float(total_reads)
		intra_bin_trans_prop = intra_bin_trans/float(total_reads)
		intra_binless_trans_prop= intra_binless_trans/float(total_reads)
		inter_bin_trans_prop= inter_bin_trans/float(total_reads)
		inter_binless_trans_prop= inter_binless_trans/float(total_reads)
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n'.format(sample, mobility, total_reads, total_trans, intra_bin_cis, intra_binless_cis, intra_bin_trans, intra_binless_trans, inter_bin_trans, inter_binless_trans, intra_bin_cis_prop, intra_binless_cis_prop, intra_bin_trans_prop, intra_binless_trans_prop, inter_bin_trans_prop, inter_binless_trans_prop))





