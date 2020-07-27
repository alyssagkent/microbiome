#normalize_trans.py

import argparse
from argparse import RawDescriptionHelpFormatter
from Bio import SeqIO
#def getOptions():
description="""normalization script, input trans file, combo file, output: """

parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-r','--rpkm', dest="rpkmhandle", action='store', required=True,  help='RPKM file for metagenomes [REQUIRED]', metavar="INFILE")
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED', metavar="OUTFILE")
parser.add_argument('-t','--trans', dest="transhandle", action='store', required=True,  help='Transfile ncol format, noeuks, reconnected? [REQUIRED]', metavar="TRANSFILE")
parser.add_argument('-s','--scf', dest="scaffold", action='store', required=True,  help='Scaffolds/Contigs [REQUIRED]', metavar="SCF")
args = parser.parse_args()
#

#get the length of the contigs
length_dict = {}
rf_dict = {}
for seq_record in SeqIO.parse(args.scaffold, "fasta"):
    seqid = seq_record.id
    l = len(seq_record)
    r = seq_record.seq.count("GATC")
    length_dict[seqid] = l
    rf_dict[seqid] = r + 1 #plus one because you have one extra fragment than the cutsites

#get rpkm information
rpkm_dict = {}
with open(args.rpkmhandle) as rpkmfile:
	header = rpkmfile.readline()
	for line in rpkmfile:
		contig = line.split(',')[0]
		rpkm = float(line.split(',')[1])
		covg = float(line.strip().split(',')[2])
		#if covg >= 80:
		rpkm = float(line.split(',')[1])
		#else:
		#	rpkm = 0
		rpkm_dict[contig] = rpkm


trans_dict = {}
with open(args.transhandle) as trans:
	for line in trans:
		ref1 = line.split('\t')[0]
		ref2 = line.split('\t')[1]
		count = float(line.strip().split('\t')[2])
		if count>=5:
			#I want a dictionary of non-double-counted reads
			try:
				#try the first combination
				trans_dict[ref1+'@'+ref2]+=count
			except KeyError:
				try:
					#if first doesn't work, try the second combination
					trans_dict[ref1+'@'+ref2]+=count
				except KeyError:
					#if the second doesn't work either then the first will be the one we use
					trans_dict[ref1+'@'+ref2]=count


scaling_factor = 10000000 #scaling factor of ten million
#go through and normalize to contig length, contig abundances and scaled by X
with open(args.outhandle,'w') as outfile:
	for key in trans_dict.keys():
		ref1 = key.split('@')[0]
		ref2 = key.split('@')[1]
		l1 = length_dict[ref1]
		l2 = length_dict[ref2]
		r1 = rpkm_dict[ref1]
		r2 = rpkm_dict[ref2]
		rf1 = rf_dict[ref1]
		rf2 = rf_dict[ref2]
		count=trans_dict[key]
		normalized = scaling_factor * count/(l1*l2*r1*r2)
		normalized_rf = scaling_factor * count/(rf1*rf2*r1*r2)
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(ref1, ref2, normalized, normalized_rf, count,l1,l2,r1,r2, rf1,rf2))
	

	
