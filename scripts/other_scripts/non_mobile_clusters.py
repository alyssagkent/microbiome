#non_mobile_clusters.py
#goal: input a fasta file of contigs, run through and spit out a fasta file without mobile contigs
#USAGE: python non_mobile_clusters.py -i Cluster_0.fasta -o Cluster_0_non_mobile.fasta -c /workdir/users/agk85/press2/

import argparse
from argparse import RawDescriptionHelpFormatter
from Bio import SeqIO
import glob

def getOptions():
	"""Get arguments"""
	description="""This script can be used to get arg versus organism capturing orgs on the 
	contig and contigs up to N links away via Hi-C reads"""
	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',dest='inhandle',action='store',required=True,type=str, help='Fasta file', metavar='INFASTA')
	parser.add_argument('-o','--output',dest='outhandle',action='store',required=True,type=str, help='Fasta file without mobile contigs', metavar='OUTFASTA')
	parser.add_argument('-c','--combo_tables',dest='combo_tables',action='store',required=True,type=str, help='Path to combo_tables', metavar='Combo_tables')
	args = parser.parse_args()
	return(args)
	
args = getOptions()
inhandle = args.inhandle
outhandle = args.outhandle
combo_tables = args.combo_tables

#define mobile
tblpaths = glob.glob(combo_tables + '/metagenomes/*_master_scf_table.txt')
is_mobile = {}
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			mge = line.split('\t')[-2]
			if mge == '.':
				mobile = 0
			if mge == 'mge':
				mobile = 1
			is_mobile[line.split('\t')[0]] = mobile

#run through and check mobility, add to rec list if it isn't mobile
recs = []
for rec in SeqIO.parse(inhandle, "fasta"):
	seqid = rec.id
	try:
		mobile = is_mobile[seqid]
		if not mobile:
			recs.append(rec)
	except KeyError:
		print(seqid, " is eukaryotic?")

#output non-mobile fastas to outhandle
SeqIO.write(recs, outhandle, 'fasta')	
