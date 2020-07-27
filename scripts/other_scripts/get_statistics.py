#get N50
#USAGE: python ~/agk/CDC/scripts/get_statistics.py genome
import sys
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import numpy
import os
import glob

type = sys.argv[1]
#paths = glob.glob('/home/agk85/agk/CDC/idba/r2_precorrection/*/scaffold_min_500.fa')
paths = glob.glob('/workdir/users/agk85/CDC2/metaspades/' + type + 's/scaffolds/*scaffold.fasta')

#file length
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


outfile = open('/workdir/users/agk85/CDC2/metaspades/metagenomes/' + type + '_statistics.txt','w')
outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format('Name','reads','scfs','total','mean','median','max','n50'))
for path in paths:
	handle = path
	lengths=[]
	for seq_record in SeqIO.parse(handle, "fasta"):
		seqlen = len(seq_record)
		lengths.append(seqlen)
	#initial stats
	scf_count = len(lengths)
	total = sum(lengths)
	mean = numpy.mean(lengths)
	median= numpy.median(lengths)
	maxvalue = max(lengths)
	half = sum(lengths)/2
	count = 0
	lengths.sort(reverse=True)
	for item in lengths:
		count = count + item
		if count > half:
			n50 = item
			break
	
	root = path.split('/')[8].split('_')[0]
	#if type == 'genome':
	#	reads = str(file_len('/workdir/data/CDC/genomes/' + root + '.derep_1.adapter.fastq')/4)
	#if type == 'metagenome':
	#	reads = str(file_len('/workdir/data/CDC/metagenomes/' + root + '.1.fastq')/4)
	outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(root,'NA',str(scf_count),str(total),str(mean),str(median),str(maxvalue),str(n50)))


outfile.close()
