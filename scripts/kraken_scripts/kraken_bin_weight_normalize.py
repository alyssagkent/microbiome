#weight normalized kraken output for contig bins
# Author: Alyssa Kent, 2019
import argparse
from argparse import RawDescriptionHelpFormatter
from copy import deepcopy
from Bio import SeqIO

#def getOptions():
description="""Program turns kraken taxonomy output """
#example python kraken_bin_weight_normalize.py -i /workdir/users/agk85/CDC2/das/B314-1/kraken/B314-1_metabat_7.contigs.fa.kraken.taxonomy.txt -f /workdir/users/agk85/CDC2/prodigal_excise/metagenomes/B314-1/B314-1_scaffold.fasta -o /workdir/users/agk85/CDC2/das/B314-1/kraken/B314-1_metabat_7.contigs.fa.kraken.weighted_besttaxid.txt -d /workdir/users/agk85/CDC2/das/B314-1/B314-1_DASTool_scaffolds2bin.txt -b B314-1_metabat_7
parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i','--in', dest="inhandle", action='store', required=True, help='Inhandle [REQUIRED]')
parser.add_argument('-f','--fasta', dest="fastahandle", action='store', required=True, help='Fastahandle [REQUIRED]')
parser.add_argument('-b','--bin', dest="binname", action='store', required=True, help='Bin_name [REQUIRED]')
parser.add_argument('-d','--das', dest="das", action='store', required=True, help='DAS contig_bin linker [REQUIRED]')
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED]')
#parser.add_argument('-g','--gtdb', dest="gtdb_clusters", action='store', required=True, help='GTDB cluster taxonomies [REQUIRED]', metavar="GTDB")

args = parser.parse_args()

length_dict={}
for rec in SeqIO.parse(args.fastahandle, 'fasta'):
	l = len(rec)
	contig = rec.id
	length_dict[contig] = l



#so you want to pick up the non-taxonomies too 
#so get a list of all of the contigs in this sample and as you add them remove them
#then go through the end and add in more length so it'll be the true percentage

#open das file that has all the associations, make a dictionary of bin_dict[bin] = list of contigs
bin_dict = {}
with open(args.das) as binfile:
	for line in binfile:
		contig = line.split('\t')[0]
		bin = line.strip().split('\t')[1]
		try:
			bin_dict[bin].append(contig)
		except KeyError:
			bin_dict[bin] = [contig]


bin = args.binname
leftover_contigs = deepcopy(bin_dict[bin])


levels = ['kingdom','phylum','class','order','family','genus','species']
level_dict = {}
for i in range(len(levels)):
	level_dict[i] = {}


total_length = 0
with open(args.inhandle) as infile:
	for line in infile:
		contig = line.split('\t')[0]
		leftover_contigs.remove(contig)
		contig_length = length_dict[contig]
		total_length += contig_length
		taxonomy = line.strip().split('\t')[1]
		taxa = taxonomy.split('|')
		for i in range(len(levels)):
			try:
				#does it exist?
				t = taxa[i]
				taxon = ';'.join(taxa[0:i+1])
				try:
					level_dict[i][taxon]+= contig_length
				except KeyError:
					level_dict[i][taxon] = contig_length
			except IndexError:
				#it doesn't have that level
				#so add the length to empty
				try:
					level_dict[i]['NA']+= contig_length
				except KeyError:
					level_dict[i]['NA'] = contig_length
	#then go through the contigs that don't exist
	for contig in leftover_contigs:
		contig_length = length_dict[contig]
		total_length += contig_length
		for i in range(len(levels)):
			#it doesn't have that level
			#so add the length to empty
			try:
				level_dict[i]['NA']+= contig_length
			except KeyError:
				level_dict[i]['NA'] = contig_length



biggest_level=-1
biggest_taxa ='NA'
biggest_proportion=0
for i in range(len(levels)):
	for taxa in level_dict[i].keys():
		if taxa != 'NA':
			prop = level_dict[i][taxa]/float(total_length)
			if prop > 0.50: #because it has to be more than 50%, then nothing else should be able to usurp it in that particular level
				if i > biggest_level:
					biggest_level = i
					biggest_taxa = taxa
					biggest_proportion = prop




contamination = 0
if biggest_level != -1:
	for taxa in level_dict[biggest_level].keys():
			if (taxa != 'NA' and taxa != biggest_taxa):
					contamination+= level_dict[biggest_level][taxa] 

contamination_prop = contamination/float(total_length)



levelsplits = ['d__','p__','c__','o__','f__','g__','s__']
if biggest_level == -1:
	biggest_taxa = ';'.join(levelsplits) + ';'


taxonomy = ['d__','p__','c__','o__','f__','g__','s__']
biggest_taxa_units = biggest_taxa.split(';')
for i in range(len(levels)):
	levelsplit = levelsplits[i]
	try:
		taxa_of_interest = biggest_taxa.split(levelsplit)[1].split(';')[0]
	except IndexError:
		taxa_of_interest = ''
	taxonomy[i] = levelsplit + taxa_of_interest
	

if biggest_level != -1:
	taxa_levelname = levels[biggest_level]
else:
	taxa_levelname = 'Null'

joined_taxonomy = '; '.join(taxonomy) + ';'
	

full_taxonomy = joined_taxonomy.replace('d__','k__')
with open(args.outhandle,'w') as outfile:
	outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(bin+'.contigs.fa.report.txt.besttaxid.contamination',bin+'.contigs',  str(100*biggest_proportion),'NA', taxa_levelname, str(biggest_level+1), str(100*contamination_prop), full_taxonomy))

			
	

