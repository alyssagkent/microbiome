#refined ARG vs. Organism traversal
# #use the best organism function from CDC one
#python arg_org_hic_refined.py -b /workdir/users/agk/gates -a 95 -p 0 -d 2
#goal of latest update
#instead of using just the best organism given by the combo_tables, instead use the best organism given by the clusters from a given cluster file
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
from best_org_nophage import best_org_nophage
from best_org import best_org

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
 	parser.add_argument('-d', '--depth', dest='depth', action='store', required=True, type=int, help='HiC depth [Required]', metavar='DEPTH') 
	parser.add_argument('-n', '--name',dest='outfilename',action='store',required=True, type=str,help='outputname B314')
	parser.add_argument('-s','--samples', nargs='+', dest='samples',help='List of patients: B314-1 B314-2 [Required]', required=True, metavar='SAMP')
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

#best_org_dicts
depth = args.depth
print(type(depth))
best_org_dict= best_org(tbldict,0)
best_org_dict_phaged = best_org_nophage(tbldict,depth)

g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping/trans_primary_ncol_' + '98' + '_withexcise_noeuks.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
#only keep the links greater than 2
g.es.select(weight_lt=2).delete()

samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

samples = args.samples

inhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/args_95_metagenomes_genomes_95.clstr'
protclusterdict = {}
fullclusterprotdict = {}
fullclusternum_map = {}
fullgenomeclusterprotdict = {}
with open(inhandle) as infile:
	for line in infile:
		if line[0] == '>':
			cluster = line.strip().split(' ')[1]
		else:
			prot =  line.split('>')[1].split('...')[0]
			samp = prot[0:6]
			protclusterdict[prot] = cluster
			if samp in samples:
				if 'NODE' not in prot:
					try:
						fullclusterprotdict[cluster].append(prot)
					except:
						fullclusterprotdict[cluster] = [prot]
				else:
					try:
						fullgenomeclusterprotdict[cluster].append(prot)
					except:	
						fullgenomeclusterprotdict[cluster] = [prot]
			if '*' in line:	
				fullclusternum_map[cluster] = prot

#i feel like this is unnecessary
clusterprotdict = {}
clusternum_map = {}
for cluster in fullclusterprotdict.keys():
	clusterprotdict[cluster] = fullclusterprotdict[cluster]
	clusternum_map[cluster] = fullclusternum_map[cluster]


#initialize argdict
mgedict = collections.defaultdict(dict)
for key in clusterprotdict.keys():
	mgedict[key] = collections.defaultdict(dict)
	for samp in samples:
		mgedict[key] = [] #list for organisms

#this one uses 
for key in clusterprotdict.keys():
	genes = clusterprotdict[key]
	for gene in genes:
		gene_samp = gene.split('_')[0]
		if gene_samp in samples:
			scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
			try:
				arg_vids = g.neighborhood(scf, order=depth, mode=ALL)
				nodes_of_interest = g.vs[arg_vids]['name']
			except ValueError:
				nodes_of_interest = [scf]
			#remove the scf (self)
			nodes_of_interest = [x for x in nodes_of_interest if x != scf]
			if 'phage' in scf:
				parentscf = gene.split('_')[0] + '_scaffold_' + gene.split('_')[1].split('|')[0]
				nodes_of_interest = [x for x in nodes_of_interest if x != parentscf] 
			for node in nodes_of_interest:
				try:
					besttaxa = best_org_dict[node]
					for besttaxon in besttaxa:
						if besttaxon != '.':
							mgedict[key].append(besttaxon)
				except KeyError:
					besttaxa = ''


#initialize argdict
basedict = collections.defaultdict(dict)
for key in clusterprotdict.keys():
	basedict[key] = []

for key in clusterprotdict.keys():
	genes = clusterprotdict[key]
	for gene in genes:      
		gene_samp = gene.split('_')[0]
		if gene_samp in samples:
			scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
			nodes_of_interest = [scf]
			for node in nodes_of_interest:
				try:
					besttaxa = best_org_dict_phaged[node]
					for besttaxon in besttaxa:
						if besttaxon != '.':
							basedict[key].append(besttaxon)
				except KeyError:
					besttaxa = ''


namedict = {}
mechdict = {}
broadmechdict = {}
deschandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/cliques/arg_v_samp_95_95_names_mech.txt'
with open(deschandle) as descfile:
	header = descfile.readline()
	for line in descfile:
		cluster = line.split('\t')[0]
		descriptions = line.split('\t')[2]
		mechanism = line.split('\t')[5]
		sub_mechanism=line.strip().split('\t')[6]
		namedict[cluster] = descriptions
		mechdict[cluster] = mechanism
		broadmechdict[cluster] = sub_mechanism


#######################separate org and arg and connection
def get_code(mgedict, key, org):
	connection = '0'
	if org in mgedict[key]:
		connection = '1'
	return connection

print("outputting combinations")

#protdict = SeqIO.to_dict(SeqIO.parse("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/arg_prot.fasta", "fasta"))
count = 0
protein_list = []
outhandle = '{0}/arg_v_org/metagenomes3/isolate_comparisons/arg_base_hic_withtaxa_{1}_{2}_{3}_{4}_{5}.tbl'.format('/workdir/users/agk85/CDC' , '95' , '98' ,str(depth), str(2),args.outfilename)
with open(outhandle,'w') as outfile:
	outfile.write("Cluster\tGenes\tGenome_genes\tNumber_genes_in_cluster\tName\tBroad_mechanism\tNumber_base_taxa\tNumber_hic_taxa\tBase_taxonomies\tHic_taxonomies\n")
	for key in clusterprotdict.keys():
		name = namedict[key]
		broadmechanism = broadmechdict[key]
		basetaxonomy_list=basedict[key]
		hictaxonomy_list = mgedict[key]
		numbasetaxa = len(set(basetaxonomy_list))
		numhictaxa = len(set(hictaxonomy_list))
		hictaxonomies =  ','.join(list(set(hictaxonomy_list)))
		if hictaxonomies == '':
			hictaxonomies = 'NA'
		basetaxonomies = ','.join(list(set(basetaxonomy_list)))
		if basetaxonomies == '':
			basetaxonomies = 'NA'
		gene = clusternum_map[key]
		genes = ','.join(clusterprotdict[key])
		try:
			genomegenes = ','.join(fullgenomeclusterprotdict[key])
		except KeyError:
			genomegenes = 'NA'
		if genomegenes !='NA':
			numgenes = str(len(clusterprotdict[key]))
			outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(key, genes,genomegenes, numgenes, name,broadmechanism, numbasetaxa,numhictaxa, basetaxonomies,hictaxonomies))


