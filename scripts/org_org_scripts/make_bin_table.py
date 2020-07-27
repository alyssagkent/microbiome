import numpy as np
import glob
import collections
from Bio import SeqIO
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script will make a list of all of the bins that are not bad bins >10% and give the weighted kraken taxonomy
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-k', '--kraken_bins', dest='kraken_bins', action='store', required=True, help='Kraken_bins [Required]')
	parser.add_argument('-c','--checkm', dest="checkm", action='store', required=True, help='checkm_info [REQUIRED]', metavar="CHECKM")
	parser.add_argument('-o','--out',dest="outhandle",action='store',required=True,help='outfile for bin table', metavar="OUTFILE")
	args = parser.parse_args()
 	return(args)

args = getOptions()
####################################
#make a table that links the bins with good quality and taxonomies
####BIN STUFF######
print('Checkm quality')
qualitydict = {}
binids = []
with open(args.checkm) as checkm:
	for line in checkm:
		binid = line.split('\t')[0].split('.contigs')[0]
		binids.append(binid)
		completion = float(line.split('\t')[2])
		contamination = float(line.split('\t')[3])
		quality = 'BAD'
		if (completion > 90 and contamination < 5):
			quality = 'HQ'
		elif (completion >= 50 and contamination < 10):
			quality = 'MQ'
		elif (completion < 50 and contamination < 10):
			quality = 'LQ'	
		qualitydict[binid] = quality

taxadict = {}
with open(args.kraken_bins) as krakenfile:
	for line in krakenfile:
		binid = line.split('.contigs.fa.report.txt.besttaxid')[0]
		taxonomy = line.strip().split('\t')[7]
		taxadict[binid]=taxonomy

with open(args.outhandle,'w') as outfile:
	for binid in binids:
		quality = qualitydict[binid]
		if quality != 'BAD': 
			sample = binid.split('_')[0].split('.')[0]
			patient = sample.split('-')[0]
			try:
				taxonomy = taxadict[binid]
			except KeyError:
				taxonomy = '.'
			if taxonomy == 'k__; p__; c__; o__; f__; g__; s__;':
				taxonomy = '.'
			outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(binid,patient,sample,quality,taxonomy))
