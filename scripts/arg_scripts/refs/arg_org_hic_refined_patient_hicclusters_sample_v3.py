#refined ARG vs. Organism traversal
# #use the best organism function from CDC one
#python arg_org_hic_refined.py -b /workdir/users/agk/gates -a 95 -p 0 -d 2
#goal of latest update
#instead of using just the best organism given by the combo_tables, instead use the best organism given by the clusters from a given cluster file
#
import numpy as np
import matplotlib.pyplot as plt
import glob
import collections
import sys
from Bio import SeqIO
import igraph
from igraph import *
import argparse
from argparse import RawDescriptionHelpFormatter
from best_org import best_org

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
 	parser.add_argument('-d', '--depth', dest='depth', action='store', required=True, type=int, help='HiC depth [Required]', metavar='DEPTH') 
 	args = parser.parse_args()
 	return(args)
 
args = getOptions()
#tbldict
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC' + '/combo_tables/metagenomes4/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()


#best_org_dict
depth = args.depth
best_org_dict = best_org(tbldict, depth)


top_args = []
with open('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/COI_clusters.txt') as coi:
	for line in coi:
		top_args.append(line.split('\t')[0])


#map the names back #this is specific to 95%
arg_name_dict = {}
cardres = {}
prot_arg_name = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95_names.txt'
with open(prot_arg_name) as f:
	for line in f:
		r = line.split('\t')[0].strip()
		if r != 'ORF_ID':
			arg_name_dict[r] = line.strip().split('\t')[1]
			cardres[r] = line.strip().split('\t')[3] + ":" + line.strip().split('\t')[2]

#make a shortened list of argnames
argnames = list(set(arg_name_dict.values()))
argnames.sort()

g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping/trans_primary_ncol_' + '98' + '_withexcise_noeuks.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
g.es.select(weight_lt=2).delete()

samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

inhandle = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes3/args_' + '95' + '_nr.fna.clstr'
protclusterdict = {}
clusterprotdict = {}
clusterpatientdict = {}
clusternum_map = {}
with open(inhandle) as infile:
	for line in infile:
		if line[0] == '>':
			cluster = line.strip().split(' ')[1]
		else:
			prot =  line.split('>')[1].split('...')[0]
			samp = prot.split('_')[0]
			pat = samp.split('-')[0]
			protclusterdict[prot] = cluster
			if samp in samples:
				try:
					clusterprotdict[cluster].append(prot)
					clusterpatientdict[cluster].append(pat)
				except:
					clusterprotdict[cluster] = [prot]
					clusterpatientdict[cluster] = [pat]
			if '*' in line:	
				clusternum_map[cluster] = prot


	

#map the names back #this is specific to 95%
arg_name_dict = {}
cardres = {}
prot_arg_name = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95_names.txt'
with open(prot_arg_name) as f:
	header = f.readline()
	for line in f:
		r = line.split('\t')[0].strip()
		if r != 'ORF_ID':
			arg_name_dict[r] = line.strip().split('\t')[1]
			try:
				cardres[protclusterdict[r]] = line.strip().split('\t')[3] + ":" + line.strip().split('\t')[2]
			except KeyError:
				continue
				print(r)
				#this means that the protein is not 


#make a shortened list of argnames
argnames = list(set(arg_name_dict.values()))
argnames.sort()


#initialize argdict
argdict = collections.defaultdict(dict)
for key in clusterprotdict.keys():
	argdict[key] = collections.defaultdict(dict)
	for samp in samples:
		argdict[key][samp] = [] #list for organisms

argsampdict = {}
for key in clusterprotdict.keys():
	argsampdict[key] = []
	genes = clusterprotdict[key]
	for gene in genes:
		gene_samp = gene.split('_')[0]
		argsampdict[key].append(gene_samp)

#argdict['B314-1_scaffold_1_1']['B316-1']['taxa1','taxa2','taxa3']

samporgs = {}
for samp in samples:
	samporgs[samp] = []

organisms = []
for key in clusterprotdict.keys():
	genes = clusterprotdict[key]
	for gene in genes:
		gene_samp = gene.split('_')[0]
		if gene_samp in samples:
			scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
			try:
				arg_vids = g.neighborhood(scf, order=int(depth), mode=ALL)
				nodes_of_interest = g.vs[arg_vids]['name']
			except ValueError:
				#this is when the scaffold isn't in the network
				nodes_of_interest = [scf]
			for node in nodes_of_interest:
				try:
					besttaxa = best_org_dict[node]
					for besttaxon in besttaxa:
						if besttaxon != '.':
							argdict[key][gene_samp].append(besttaxon)
							organisms.append(besttaxon)
							samporgs[gene_samp].append(besttaxon)
				except KeyError:
					besttaxa = ''

sampleorgs = {}
for samp in samporgs:
	sampleorgs[samp] = list(set(samporgs[samp]))


#do i still need this?
#get all of the organisms...make a set
orglist = list(set(organisms))
orglist.sort()


##############################################
print("onto org arg counting")
delims = ['s__','g__', 'f__','o__','c__','p__','k__']
delim_names=['species','genus','family','order','class','phylum','kingdom']

delims = ['s__']
delim_names = ['species']

patients = [samp.split('-')[0] for samp in samples]
patients = list(set(patients))
print(patients)
for i in range(len(delims)):
	arghisthandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/arg_orgcounts_allsamples_' + delim_names[i] + str(args.depth) + '.txt'
	count = 0
	delim = delims[i]
	with open(arghisthandle, 'w') as argout:
		for arg in clusterprotdict.keys():
			for samp in samples:
				count = count + 1
				orgs = argdict[arg][samp]
				level_orgs = list(set([org.split(delim)[1].split(';')[0] for org in orgs]))
				try:
					level_orgs.remove('')
				except:
					a = 1
				l = len(level_orgs)
				argname = arg_name_dict[clusternum_map[arg]]
				orgstring = ','.join(level_orgs)
				patient = samp.split('-')[0]
				tp = samp.split('-')[1]
				argout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(str(count), samp, patient, tp, arg, argname, str(l), orgstring))


