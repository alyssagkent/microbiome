#contig_vs_bin.py
import argparse
from argparse import RawDescriptionHelpFormatter
from copy import deepcopy

#def getOptions():
description="""Program creates a contig vs. bin table summing all contacts between the contig and the bin of interest, not counting contacts to itself (Only trans)"""

parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="INFILE")
parser.add_argument('-l','--hic', dest="hichandle", action='store', required=True, help='HI-C file [REQUIRED]', metavar="HIC_FILE")
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED]', metavar="OUTFILE")
parser.add_argument('-t','--ct', dest="combotable", action='store', required=True, help='Combotable [REQUIRED]', metavar="OUTFILE")


args = parser.parse_args()

sample = args.combotable.split('/')[-1].split('_master_scf_table.txt')[0]
bindict = {}
with open(args.binhandle) as binfile:
	for line in binfile:
		contig = line.split('\t')[0]
		binid = line.strip().split('\t')[1]
		bindict[contig] = binid

contigs = []
with open(args.combotable) as infile:
	header = infile.readline()
	for line in infile:
		contig = line.split('\t')[0]
		contigs.append(contig)


#so want a few things
num_hic_trans = 0
num_hic_supportive = 0
num_hic_novel = 0
num_hic_inclusive = 0
num_hic_lonelyislands = 0
withinbin_supported_contigs = []
betweenbin_supported_contigs = []
lonelybin_supported_contigs =[]
lonelyisland_supported_contigs = []

with open(args.hichandle) as hicfile:
	for line in hicfile:
		#Get the contigs linking
		contig1 = line.split('\t')[1]
		contig2 = line.split('\t')[4]
		if contig1 != contig2:
			num_hic_trans += 1
			try:
				bin1 = bindict[contig1]
				try:
					bin2 = bindict[contig2]
					if bin1 == bin2:
						num_hic_supportive +=1
						withinbin_supported_contigs.append(contig1)
						withinbin_supported_contigs.append(contig2)
					else:
						#they are in different bins
						num_hic_novel +=1
						betweenbin_supported_contigs.append(contig1)
						betweenbin_supported_contigs.append(contig2)
				except KeyError:
					#bin1 but not bin2
					num_hic_inclusive += 1
					lonelybin_supported_contigs.append(contig1)
			except KeyError:
				try:
					bin2 = bindict[contig2]
					#bin2 but not bin1
					num_hic_inclusive += 1
					lonelybin_supported_contigs.append(contig2)
				except KeyError:
					#neither annotated
					num_hic_lonelyislands +=1
					lonelyisland_supported_contigs.append(contig1)
					lonelyisland_supported_contigs.append(contig2)



a = str(len(set(withinbin_supported_contigs)))
b = str(len(set(betweenbin_supported_contigs)))
c = str(len(set(lonelybin_supported_contigs)))
d = str(len(set(lonelyisland_supported_contigs)))

a1 = str(num_hic_trans)
a2 = str(num_hic_supportive)
a3 = str(num_hic_novel)
a4 = str(num_hic_inclusive)
a5 = str(num_hic_lonelyislands)


with open(args.outhandle,'w') as outfile:
	header = 'Sample\tNum_hicreads_trans\tNum_hicreads_supportive\tNum_hicreads_novel\tNum_hicreads_inclusive\tNum_hicreads_lonelyislands\tNum_contigs_withinbin_supported\tNum_contigs_betweenbin_supported\tNum_contigs_lonelybin_supported\tNum_contigs_lonelyisland_supported\tNum_contigs_unsupported\n'
	outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(sample, a1,a2,a3,a4,a5,a,b,c,d))
