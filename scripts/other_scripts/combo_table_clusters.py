#combo_table_clusters.py
#USAGE:python ~/agk/CDC2/scripts/combo_table_clusters.py -i /workdir/users/agk85/press2/combo_tables/metagenomes/ProxiMeta-1_master_scf_table.txt -o /workdir/users/agk85/press2/combo_tables/metagenomes/ProxiMeta-1_master_scf_cluster_table.txt -d /workdir/users/agk85/press2/das/ProxiMeta-1/ProxiMeta-1_DASTool_bins -n /workdir/users/agk85/press2/das/ProxiMeta-1/ProxiMeta-1_DASTool_non_mobile_bins
import argparse
from argparse import RawDescriptionHelpFormatter
from Bio import SeqIO
import glob

def getOptions():
	"""Get arguments"""
	description="""This script will add das cluster names, hierarchical cluster names, non_mobile_DAS cluster names to combo_table """
	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',dest='inhandle',action='store',required=True,type=str, help='Combo_table', metavar='INFILE')
	parser.add_argument('-o','--output',dest='outhandle',action='store',required=True,type=str, help='Fasta file without mobile contigs', metavar='OUTFILE')
	parser.add_argument('-d','--das_folder_file',dest='das',action='store',required=True,type=str, help='DAS_folder', metavar='DAS')
	parser.add_argument('-n','--das_non_mobile_folder_file',dest='das_non_mobile',action='store',required=True,type=str, help='DAS_non_mobile_folder', metavar='DAS_NON_MOBILE')
	args = parser.parse_args()
	return(args)
	
args = getOptions()
inhandle = args.inhandle
outhandle = args.outhandle
das = args.das
das_non_mobile = args.das_non_mobile

das_dict = {}
#for each file in the given folder, assign the file to the contigs it has
for fasta in glob.glob(das + '/*.fa'):
	base = fasta.split('/')[-1].split('.fa')[0]
	with open(fasta) as infile:
		for line in infile:
			if line[0] == '>':
				contig = line.strip().split('>')[1]
				das_dict[contig] = base


das_non_mobile_dict = {}
#for each file in the given folder, assign the file to the contigs it has
for fasta in glob.glob(das_non_mobile + '/*.fa'):
	base = fasta.split('/')[-1].split('.fa')[0]
	with open(fasta) as infile:
		for line in infile:
			if line[0] == '>':
				contig = line.strip().split('>')[1]
				das_non_mobile_dict[contig] = base



print(das_non_mobile_dict)

with open(inhandle) as infile, open(outhandle,'w') as outfile:
	header = infile.readline()
	newheader = header.strip() + '\tDAS\tDAS_non_mobile\n'
	outfile.write(newheader)
	for line in infile:
		contig = line.split('\t')[0]
		try:
			das_cluster = das_dict[contig]
		except KeyError:
			das_cluster = '.'
		try:
			das_non_mobile_cluster = das_non_mobile_dict[contig]
		except KeyError:
			das_non_mobile_cluster = '.'
		outfile.write(line.strip() + '\t' + das_cluster + '\t' + das_non_mobile_cluster + '\n')

