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

best_org_dict= best_org(tbldict,depth)

g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping/trans_primary_ncol_' + '98' + '_withexcise_noeuks.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
#only keep the links greater than 2
g.es.select(weight_lt=2).delete()

samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)



#only keep the clusters that are machinery
machinerydict = {}
machinery_handle = '/workdir/users/agk85/CDC/tables/metagenomes3/mge_95_nr.fna.clstr.desc'
with open(machinery_handle) as mach:
	header = mach.readline()
	for line in mach:
		cluster = line.split('\t')[0]
		phage_mach = line.split('\t')[3]
		machinerydict[cluster]=phage_mach

inhandle = '/workdir/users/agk85/CDC' + '/tables/metagenomes3/mge_' + '95' + '_nr.fna.clstr'
protclusterdict = {}
fullclusterprotdict = {}
fullclusternum_map = {}
with open(inhandle) as infile:
	for line in infile:
		if line[0] == '>':
			cluster = line.strip().split(' ')[1]
		else:
			prot =  line.split('>')[1].split('...')[0]
			samp = prot.split('_')[0]
			protclusterdict[prot] = cluster
			if samp in samples:
				try:
					fullclusterprotdict[cluster].append(prot)
				except:
					fullclusterprotdict[cluster] = [prot]
			if '*' in line:	
				fullclusternum_map[cluster] = prot

clusterprotdict = {}
clusternum_map = {}
for cluster in fullclusterprotdict.keys():
	phagemach = machinerydict[cluster]
	if phagemach == '1':
		clusterprotdict[cluster] = fullclusterprotdict[cluster]
		clusternum_map[cluster] = fullclusternum_map[cluster]


#initialize argdict
mgedict = collections.defaultdict(dict)
for key in clusterprotdict.keys():
	mgedict[key] = collections.defaultdict(dict)
	for samp in samples:
		mgedict[key] = [] #list for organisms

mgesampdict = {}
for key in clusterprotdict.keys():
	mgesampdict[key] = []
	genes = clusterprotdict[key]
	for gene in genes:
		gene_samp = gene.split('_')[0]
		mgesampdict[key].append(gene_samp)


organisms = []
for key in clusterprotdict.keys():
	genes = clusterprotdict[key]
	for gene in genes:
		gene_samp = gene.split('_')[0]
		if gene_samp in samples:
			scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
			try:
				arg_vids = g.neighborhood(scf, order=args.depth, mode=ALL)
				nodes_of_interest = g.vs[arg_vids]['name']
			except ValueError:
				nodes_of_interest = [scf]
			for node in nodes_of_interest:
				try:
					besttaxa = best_org_dict[node]
					for besttaxon in besttaxa:
						if besttaxon != '.':
							mgedict[key].append(besttaxon)
							organisms.append(besttaxon)
				except KeyError:
					besttaxa = ''

#get all of the organisms...make a set
orglist = list(set(organisms))
orglist.sort()

#mge machinery
descdict = {}
deschandle = '/workdir/users/agk85/CDC/tables/metagenomes3/mge_95_nr.fna.clstr.desc'
with open(deschandle) as descfile:
	header = descfile.readline()
	for line in descfile:
		cluster = line.split('\t')[0]
		descriptions = line.split('\t')[1]
		descdict[cluster] = descriptions


sourcedict = {}
phagesourcedict = {}
sourcehandle = '/workdir/users/agk85/CDC/tables/metagenomes3/mge_95_nr.fna.clstr.tbl.annot'
with open(sourcehandle) as sourcefile:
	header = sourcefile.readline()
	for line in sourcefile:
		cluster = line.split('\t')[0]
		source = line.split('\t')[3]
		phage = line.split('\t')[5]
		sourcedict[cluster] = source
		phagesourcedict[cluster] = phage
#######################separate org and arg and connection
def get_code(mgedict, key, org):
	connection = '0'
	if org in mgedict[key]:
		connection = '1'
	return connection

print("outputting combinations")

protdict = SeqIO.to_dict(SeqIO.parse("/workdir/users/agk85/CDC/tables/metagenomes3/mobile_genes.fna", "fasta"))

count = 0
protein_list = []
#figure out if you want to keep the lines that only have the organism...about 1 million lines
outhandle = '{0}/tables/metagenomes3/mge_machinery_hic_phage_withtaxa_{1}_{2}_{3}_{4}.tbl'.format('/workdir/users/agk85/CDC' , '95' , '98' ,str(args.depth), str(2))
with open(outhandle,'w') as outfile:
	outfile.write("Cluster\tGene\tSource\tNumber_taxa\tTaxonomies\tSpecies\tPfams\n")
	for key in clusterprotdict.keys():
		pfams = descdict[key]
		source = sourcedict[key]
		phage = phagesourcedict[key]
		taxonomy_list = mgedict[key]
		numtaxa = len(set(taxonomy_list))
		if (numtaxa >0 and int(phage)>0):
			species_list = []
			for taxonomy in taxonomy_list:
				taxa = taxonomy.split('s__')[1].split(';')[0]
				if taxa != '':
					species_list.append(taxa)
			taxonomies =  ','.join(list(set(taxonomy_list)))
			species = ','.join(list(set(species_list)))
			numspecies = len(set(species_list))
			for gene in clusterprotdict[key]:
				rec = protdict[gene]
				protein_list.append(rec)
				outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(key, gene, source, numtaxa, numspecies, taxonomies, species, pfams))

seqhandle = '{0}/tables/metagenomes3/mge_machinery_hic_phage_withtaxa_{1}_{2}_{3}_{4}.fna'.format('/workdir/users/agk85/CDC' , '95' , '98' , str(args.depth),str(2))
SeqIO.write(protein_list, seqhandle,"fasta")

