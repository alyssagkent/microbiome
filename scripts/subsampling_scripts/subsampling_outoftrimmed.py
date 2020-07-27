#rarefaction_counting.py

#the goal here is to iterate over each of the samples, iterate over all of the reads in that sample
#add the reads to the dictionary that is keeping track of the arg-org connections
#every time that the threshold 2 and 5 is passed then the connection is added to the respective lists----in both directions
#at your modulo counter then assess the number of connections that make the threshold cutoff 2 and 5 
#, then when you modulo you uniquie them and then divide by 2 because both ways 
#(do this at the same time and update each so you don't have to doubly do this)
#print out counter\tsample\tpatient\tthreshold\treads\tnum_connections

import numpy as np
import collections
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-gc','--genecluster',dest='genecluster',action='store',required=True,type=str, help='Genes clstr table [Required]', metavar='GENEFILE')
	parser.add_argument('-k', '--kraken_bins', dest='kraken_bins', action='store', required=True, help='KRAKEN_BINS [Required]')
	parser.add_argument('-c','--checkm', dest="checkm", action='store', required=True, help='checkm_info [REQUIRED]', metavar="CHECKM")
	parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="BINFILE")
	parser.add_argument('-l','--hic', dest="hic", action='store', required=True,  help='hic file [REQUIRED]', metavar="HICFILE")
	parser.add_argument('-o1','--out2', dest="outhandle2", action='store', required=True,  help='outconnections2 [REQUIRED]', metavar="OUTFILE2")
	parser.add_argument('-o2','--out5', dest="outhandle5", action='store', required=True,  help='outconnections5 [REQUIRED]', metavar="OUTFILE5")
	parser.add_argument('-o3','--outcon2', dest="outconbin2", action='store', required=True,  help='outconbin2 [REQUIRED]', metavar="OUTCONBIN2")
	parser.add_argument('-o4','--outcon5', dest="outconbin5", action='store', required=True,  help='outconbin5 [REQUIRED]', metavar="OUTCONBIN5")
	parser.add_argument('-w','--window', dest="window", action='store', required=True, type=int,  help='window [REQUIRED]', metavar="WINDOW")
	parser.add_argument('-r','--reads', dest="reads", action='store', required=True,  help='reads [REQUIRED]', metavar="READS")

	args = parser.parse_args()
 	return(args)

args = getOptions()
####################################
window = args.window
print(window)
print(type(window))
####BIN STUFF######
print('Checkm quality')
qualitydict = {}
cluster_length = {}
samplelist = []
with open(args.checkm) as checkm:
	for line in checkm:
		bin = line.split('\t')[0].split('.contigs')[0]
		sample = bin.split('.')[0].split('_')[0]
		samplelist.append(sample)
		completion = float(line.split('\t')[2])
		contamination = float(line.split('\t')[3])
		length = line.strip().split('\t')[5]
		quality = 'BAD'
		if (completion > 90 and contamination < 5):
			quality = 'HQ'
		elif (completion >= 50 and contamination < 10):
			quality = 'MQ'
		elif (completion < 50 and contamination < 10):
			quality = 'LQ'	
		qualitydict[bin] = quality
		cluster_length[bin] = length

samples = list(set(samplelist))
samples.sort()

#print('Initialize orgfulldict')
#orgfulldict = {}
#orgdict = {}
#for sample in samples:
#	orgfulldict[sample] = []
#	orgdict[sample] = []


print('bindict')
#organismsfull = []
bindict = {}
with open(args.kraken_bins) as krakenfile:
	for line in krakenfile:
		bin = line.split('.contigs.fa.report.txt.besttaxid')[0]
		sample = bin.split('_')[0].split('.')[0]
		taxonomy = line.strip().split('\t')[7]
		#orgfulldict[sample].append(taxonomy)
		bindict[bin] = taxonomy
		#organismsfull.append(taxonomy)



print('bincontigdict')
bincontigdict = {}
contigbindict = {}
contigs = []
with open(args.binhandle) as binfile:
	for line in binfile:
		bin = line.strip().split('\t')[1]
		contig = line.split('\t')[0]
		binquality = qualitydict[bin]
		if binquality != 'BAD':
			contigs.append(contig)
			contigbindict[contig] = bin
			try:
				bincontigdict[bin][contig] = 0
			except KeyError:
				bincontigdict[bin] = {}
				bincontigdict[bin][contig] = 0

print(type(bincontigdict))

#make the orgdict unique
#orgs = list(set(organismsfull))
#orgs.sort()
#for sample in samples:
#	orgdict[sample] = list(set(orgfulldict[sample]))


print('Create the genedict')
genedict = {}
with open(args.genecluster) as genefile:
	for line in genefile:
		cluster = line.split('\t')[0]
		genelist = line.split('\t')[2].split(',')
		for gene in genelist:
			contig = '_'.join(gene.split('_')[0:3])
			try:
				genedict[contig].append(cluster)
			except KeyError:
				genedict[contig] = [cluster]


print('Filling connections')
connections_2 = [] #these are arg-org connections 
connections_5 = []

conbin_2 = [] # these are contig vs. bin connections
conbin_5 = []


#bindict is global
#genedict is global
def update_connections(contig, bin, conbin, connections):
	conbin.append((contig, bin)) #updating the conbin
	try:
		taxonomy = bindict[bin]
		#make sure it has a species
		species = taxonomy.split('s__')[1].split(';')[0]
		if species != '': 
			geneids = genedict[contig]
			for geneid in geneids:
				connections.append((geneid, taxonomy))
	except KeyError:
		#no taxonomy or no geneids
		a = 1	


#you have to somehow load the connections that are already present because of residency
#iterate through contigs
#get the bin of the contig
#get the taxonomy of the bin
#update
print('loading residents')
#these should only be the contigs that have good bins
for contig in contigs:
	try:
		bin = contigbindict[contig]
		update_connections(contig, bin, conbin_2, connections_2)
		update_connections(contig, bin, conbin_5, connections_5)
	except KeyError:
		a = 1
	if 'phage' in contig:
		print('yeah im phaging')
		unphaged = contig.split('_')[0] + '_scaffold_' + contig.split('|')[0].split('_')[1]
		try:
			bin1 = contigbindict[contig]
			update_connections(unphaged, bin1, conbin_2, connections_2)
			update_connections(unphaged, bin1, conbin_5, connections_5)
			print('yeah im phaging1')
		except KeyError:
			a = 1
		try:
			bin2 = contigbindict[unphaged]
			update_connections(contig, bin2, conbin_2, connections_2)
			update_connections(contig, bin2, conbin_5, connections_5)
			print('yeah im phaging2')
		except KeyError:
			a = 1
			
resident_connections=len(set(connections_2))
resident_conbin = len(set(conbin_2))

print(resident_connections)
print(len(connections_2))

#ok make a dictionary for the trans reads to try
hicreaddict = {} 
with open(args.hic) as hicfile:
	for line in hicfile:
		readid = line.split('\t')[0]
		hicreaddict[readid] = line


sample = args.kraken_bins.split('/')[-1].split('_all_kraken_weighted.txt')[0]
patient = sample.split('-')[0]
count = 0
print('starting the interation through the hic')
header = 'sample\tpatient\tthreshold\treads\tunique_connections\thic_connections\n'
with open(args.reads) as readfile, open(args.outhandle2,'w') as outfile2, open(args.outhandle5,'w') as outfile5, open(args.outconbin2,'w') as outconbin2, open(args.outconbin5,'w') as outconbin5:
	outfile2.write(header)
	outfile5.write(header)
	outconbin2.write(header)
	outconbin5.write(header)
	for line1 in readfile:
		if line1[0] == '@':
			count += 1
			readid = line1.split(' ')[0].split('@')[1]
			try:
				line = hicreaddict[readid]
				#Get the contigs linking
				if 'trans' in args.hic:
					contig1 = line.split('\t')[2]
					contig2 = line.split('\t')[8]
				else:
					contig1 = line.split('\t')[1]
					contig2 = line.split('\t')[4]
				if contig1 != contig2:
					try:
						bin1 = contigbindict[contig1] #not all of them have bins
						try:
							bincontigdict[bin1][contig2] +=1 #maybe you haven't initialized it
							if bincontigdict[bin1][contig2] == 2:
								update_connections(contig2, bin1, conbin_2, connections_2)
							if bincontigdict[bin1][contig2] == 5:
								update_connections(contig2, bin1, conbin_5, connections_5)
						except KeyError:	
							bincontigdict[bin1][contig2] = 1
					except KeyError:
						a = 1
					#same thing but in reverse #are we doubling everything here??? think about that later!!!
					try:
						bin2 = contigbindict[contig2]
						try:
							bincontigdict[bin2][contig1] +=1
							if bincontigdict[bin2][contig1] == 2:
								update_connections(contig1, bin2, conbin_2, connections_2)
							if bincontigdict[bin2][contig1] == 5:
								update_connections(contig1, bin2, conbin_5, connections_5)
						except KeyError:
							bincontigdict[bin2][contig1] = 1 # you know it hasn't passed the thresh
					except KeyError:
						a = 1
					#handle the iterations
			except KeyError:
				#not a trans hic read that passed threshold
				a = 1
			if (count % window) == 0:
				print(count)
				num_connections_2 = len(set(connections_2))
				num_connections_5 = len(set(connections_5))
				num_conbin_2 = len(set(conbin_2))
				num_conbin_5 = len(set(conbin_5))
				outfile2.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(sample, patient, '2', str(count), num_connections_2, num_connections_2-resident_connections))
				outfile5.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(sample, patient, '5', str(count), num_connections_5,num_connections_5-resident_connections))
				outconbin2.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(sample, patient, '2', str(count), num_conbin_2,num_conbin_2-resident_conbin))
				outconbin5.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(sample, patient, '5', str(count), num_conbin_5,num_conbin_5-resident_conbin))




		
