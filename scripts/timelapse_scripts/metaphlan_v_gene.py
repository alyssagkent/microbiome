import numpy as np
import glob
import collections
from Bio import SeqIO
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-gc','--genecluster',dest='genecluster',action='store',required=True,type=str, help='Genes clstr file [Required]', metavar='ARGPID')
	parser.add_argument('-g','--genetype',dest='genetype',action='store',required=True,type=str, help='Gene [Required]', metavar='GENE')
	parser.add_argument('-m', '--minreads', dest='minreads', action='store', required=True, type=int, help='Min reads [Required]')
	parser.add_argument('-c','--connections', dest="connections", action='store', required=True, help='Contig vs bin file [REQUIRED]', metavar="CONNECTIONS")
	parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="INFILE")
	args = parser.parse_args()
 	return(args)

args = getOptions()
####################################
#ok simpler you just want to know for each patient, out of all times when there is a connection, how consistent is it across the timepoints

####BIN STUFF######

samplelist = []
with open(args.binhandle) as binfile:
	for line in binfile:
		binid,patient,sample,quality,taxonomy = line.strip().split('\t')
		samplelist.append(sample)

samples = list(set(samplelist))
samples.sort()

levels = ['k__','p__','c__','o__','f__','g__', 's__']

print('Initialize connectiondict')
connection_dict = {}
for sample in samples:
	connection_dict[sample] = {}
	for level in levels:
		connection_dict[sample][level] = []

#B314    B314-1  arg     2       5329    k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia_coli;
with open(args.connections) as infile:
	header = infile.readline()
	for line in infile:
		patient,sample,genetype,minreads,geneid,taxonomy = line.strip().split('\t')
		for level in levels:
			taxon = taxonomy.split(level)[1].split(';')[0]
			if taxon != '':
				connection_dict[sample][level].append((geneid,taxon))

#get the arg rpkm
rpkmdict = {}
generepdict = {}
genenamedict = {}
if genetype == 'arg':
	with open('/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech.txt') as infile:
		for line in infile:
			cluster = line.split('\t')[0]
			repgene = line.split('\t')[1]
			name = line.split('\t')[3] #give it the card name unless its empty
			if name == 'NA':
				name = line.split('\t')[2]
			genenamedict[cluster] = name
			generepdict[repgene] = cluster
	
	with open('/workdir/users/agk85/CDC2/args/mapping/bwa_alignments_99_99/arg_v_samp_99_99.txt') as rpkmfile:
		header = rpkmfile.readline()
		geneids = header.strip().split(',')
		sid = geneids.pop(0)
		for line in rpkmfile:
			rpkms = line.strip().split(',')
			sample = rpkms.pop(0)
			rpkmdict[sample] = {}
			for repgene, rpkm in zip(geneids, rpkms):
				clusterid = generepdict[repgene]
				rpkmdict[sample][clusterid] = rpkm

#get the taxa abundance
abunddict = {}
for sample in samples:
	abunddict[sample] = {}

with open('/workdir/users/agk85/CDC2/das/all_bintables_metaphlan.txt') as infile:
	for line in infile:
		binid,patient,sample,quality,taxonomy,k,p,c,o,f,g,s,best_binid_level,best_abund_level,diff = line.strip().split('\t')
		taxabunds = [k,p,c,o,f,g,s]
		if taxonomy == '.':
			taxonomy = 'k__; p__; c__; o__; f__; g__; s__;'
		for i in range(len(levels)):
			level = levels[i]
			taxabund = taxabunds[i]
			taxon = taxonomy.split(level)[1].split(';')[0]
			abunddict[sample][taxon] = taxabund


timepoint_dict = {"B314-1":"1","B314-2":"2","B314-3":"3","B314-4":"4","B316-1":"1","B316-2":"2","B316-3":"3","B316-4":"4","B316-5":"5","B316-6":"6","B316-7":"7","B320-1":"1","B320-2":"2","B320-3":"3","B320-5":"4","B331-1":"1","B331-2":"2","B331-3":"3","B331-4":"4","B335-1":"1","B335-2":"2","B335-3":"3","B357-1":"1","B357-2":"2","B357-3":"3","B357-4":"4","B357-5":"5","B357-6":"6","B370-1":"1","B370-2":"2","B370-3":"3","B370-4":"4","B370-5":"5","US3-8":"1","US3-10":"2","US3-12":"3","US3-14":"4","US3-16":"5","US8-1":"1","US8-2":"2","US8-3":"3","US8-4":"4","US8-5":"5"}

outhandle = '/workdir/users/agk85/CDC2/bins/gene_abundances/metaphlan_v_generpkm_{0}_org_{1}.txt'.format(args.genetype, str(args.minreads))
with open(outhandle, 'w') as outfile:
	header = 'patient\tsample\ttimepoint\tgenetype\tminthresh\tlevel\tclusterid\ttaxon\tgenerpkm\ttaxon_abund\n'
	outfile.write(header)
	for level in levels:
		for sample in samples:
			for connection in connection_dict[sample][level]:
				clusterid = connection[0]
				taxon = connection[1]
				timepoint = timepoint_dict[sample]
				patient = sample.split('-')[0]
				generpkm= rpkmdict[sample][clusterid]
				taxabund = abunddict[sample][taxon]
				if taxabund != 'NA':
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(patient, sample, timepoint, args.genetype, str(args.minreads),level,clusterid, taxon, generpkm,taxabund))
