# goal given a list of sampleIDs and cluster ids return the contig +/- 1000 bp of the gene. of interest
python gene_context.py genelist.txt /workdir/users/agk85/CDC2/bins/gene_abundances/gene_contexts.txt /workdir/users/agk85/CDC2/args/ 
import sys
from Bio import SeqIO

inhandle = sys.argv[1]
outhandle = sys.argv[2]
arghandle = sys.argv[3]

with open(inhandle) as infile:
	for line in infile:
		sample = line.split('\t')[0]
		clusterid = line.strip().split('\t')[1]
		
