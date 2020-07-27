#goal keep track of the number of reads linking each...
#so we might have to handle the thing about min 2 reads somehow

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

depth = 1
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
		argdict[key][samp] = {} #list for organisms

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
				arg_vids = g.neighborhood(scf, order=depth, mode=ALL)
				nodes_of_interest = g.vs[arg_vids]['name']
				for node in nodes_of_interest:
					try:
						besttaxa = best_org_dict[node]
						for besttaxon in besttaxa:
							if besttaxon != '.':
								readcount = g[scf, node]
								try:
									argdict[key][gene_samp][besttaxon]+= readcount
								except KeyError:
									argdict[key][gene_samp][besttaxon] = readcount
								organisms.append(besttaxon)
								samporgs[gene_samp].append(besttaxon)
					except KeyError:
						besttaxa = ''
			except ValueError:
				a = 1

sampleorgs = {}
for samp in samporgs:
	sampleorgs[samp] = list(set(samporgs[samp]))

#get all of the organisms...make a set
orglist = list(set(organisms))
orglist.sort()


#######################separate org and arg and connection
header = 'Count\tCluster\tARG_name\tOrganism'
for samp in samples:
	header = header + '\t' + samp

header = header +'\n'

print("outputting combinations")

count = 0
#figure out if you want to keep the lines that only have the organism...about 1 million lines
outhandle = '{0}/arg_v_org/metagenomes3/arg_org_read_abundance_{1}_{2}_{3}_{4}.tbl'.format('/workdir/users/agk85/CDC' , '95' , '98' ,str(depth), str(2))
with open(outhandle,'w') as outfile:
	outfile.write(header)
	for key in clusterprotdict.keys():
		argname = arg_name_dict[clusternum_map[key]]
		for organism in orglist:
			count += 1
			linesum = 0
			tentativeline = ''
			tentativeline = '{0}\t{1}\t{2}\t{3}'.format(str(count), key, argname, organism)
			for samp in samples:
				try:
					readcount = argdict[key][samp][organism]
					linesum = linesum + readcount
				except KeyError:
					readcount = 0
				tentativeline = tentativeline +'\t' + str(readcount)
			if linesum > 0:	
				tentativeline = tentativeline + '\n'
				outfile.write(tentativeline)


