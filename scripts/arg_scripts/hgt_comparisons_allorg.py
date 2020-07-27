#smillie esque taxa comparisons at genus level first
#improve to looking at different taxonomic levels
import numpy as np
import glob
import collections
import sys
from Bio import SeqIO
import igraph
from igraph import *
import argparse
from argparse import RawDescriptionHelpFormatter
from best_org import best_org 
 
#tbldict
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC2' + '/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()

depth = 1
folder = 'CDC2'
argpid = '99'
hicpid = '0'
combohandle = '/workdir/users/agk85/CDC2/ComboDesign.txt'
minreads=2

best_org_dict= best_org(tbldict,depth)

#map the names back #this is specific to 95%
arg_name_dict = {}
cardres = {}
prot_arg_name = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes/mapping/bwa_alignments_99_99/arg_v_samp_99_99_names.txt'
with open(prot_arg_name) as f:
	for line in f:
		r = line.split('\t')[0].strip()
		if r != 'ORF_ID':
			arg_name_dict[r] = line.strip().split('\t')[1]
			cardres[r] = line.strip().split('\t')[3] + ":" + line.strip().split('\t')[2]

#make a shortened list of argnames
argnames = list(set(arg_name_dict.values()))
argnames.sort()

g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping/trans_primary_ncol_' + '0' + '_withexcise_noeuks.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
#only keep the links greater than 2
g.es.select(weight_lt=2).delete()


samples = []
patients = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)
	patients.append(samp.split('-')[0])

patients = list(set(patients))
patients.sort()

inhandle = '/workdir/users/agk85/CDC' + '/arg_v_org/metagenomes3/args_' + '95' + '_nr.fna.clstr'
protclusterdict = {}
clusterprotdict = {}
clusternum_map = {}
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
					clusterprotdict[cluster].append(prot)
				except:
					clusterprotdict[cluster] = [prot]
			if '*' in line:	
				clusternum_map[cluster] = prot

#initialize argdict
argdict = collections.defaultdict(dict)
for key in clusterprotdict.keys():
	argdict[key] = collections.defaultdict(dict)
	for pat in patients:
		argdict[key][pat] = [] #list for organisms

argpatdict = {}
for key in clusterprotdict.keys():
	argpatdict[key] = []
	genes = clusterprotdict[key]
	for gene in genes:
		gene_pat = gene.split('-')[0]
		argpatdict[key].append(gene_pat)

patorgs = {}
for pat in patients:
	patorgs[pat] = []

organisms = []
for key in clusterprotdict.keys():
	genes = clusterprotdict[key]
	for gene in genes:
		gene_samp = gene.split('_')[0]
		gene_patient = gene.split('-')[0]
		if gene_samp in samples:
			scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
			try:
				arg_vids = g.neighborhood(scf, order=depth, mode=ALL)
				nodes_of_interest = g.vs[arg_vids]['name']
			except ValueError:
				nodes_of_interest = [scf]
			for node in nodes_of_interest:
				try:
					besttaxa = best_org_dict[node]
					for besttaxon in besttaxa:
						if besttaxon != '.':
							argdict[key][gene_patient].append(besttaxon)
							organisms.append(besttaxon)
							#patorgs[gene_patient].append(besttaxon)
				except KeyError:
					besttaxa = ''


patientorgs = {}
for node in best_org_dict:
	besttaxa = best_org_dict[node]
	for besttaxon in besttaxa:
		if besttaxon != '.':
			pat = node.split('-')[0]
			patorgs[pat].append(besttaxon)


patientorgs = {}
for pat in patorgs:
	patientorgs[pat] = list(set(patorgs[pat]))

#get all of the organisms...make a set
orglist = list(set(organisms))
orglist.sort()

#################################################################ok so let's begin
#taxalists
delims = ['g__','f__','o__','c__','p__','k__']
delimnames = ['genus','family','order','class','phylum','kingdom']
taxalists = collections.defaultdict(dict)
for patient in patients:
	taxalists[patient] = collections.defaultdict(dict)
	for delimname in delimnames:
		taxalists[patient][delimname]= collections.defaultdict(dict)


alltaxa = collections.defaultdict(dict)

for i in range(len(delimnames)):
	delimname = delimnames[i]
	delim = delims[i]
	alltaxa[delimname] = []
	for patient in patients:
		for org in patientorgs[patient]:
			if org.split('s__')[1].split(';')[0] != '':
				taxa = org.split(delim)[1].split(';')[0]
				if taxa != '':
					alltaxa[delimname].append(taxa)
					try:
						taxalists[patient][delimname][taxa].append(org)
					except AttributeError:
						taxalists[patient][delimname][taxa] = [org]


for delimname in delimnames:
	alltaxa[delimname] = list(set(alltaxa[delimname]))
	alltaxa[delimname].sort()

#comparing two people
finished = []
outhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/hgt_comparisons/hgt_comparisons_alllevels_allorgs.txt'
with open(outhandle, 'w') as outfile:
	for patient1 in patients:
		for patient2 in patients:
			if (patient2,patient1) not in finished:
				finished.append((patient1,patient2))
				for i in range(len(delimnames)):
					delimname = delimnames[i]
					linkedcount = 0
					allcount = 0
					for taxa in alltaxa[delimname]:
						t1 = taxa in taxalists[patient1][delimname].keys()
						t2 = taxa in taxalists[patient2][delimname].keys()
						if (t1 and t2):
							for species1 in taxalists[patient1][delimname][taxa]:
								for species2 in taxalists[patient2][delimname][taxa]:
									if i>0:
										#assess whether or not they are the same i-1 taxa 
										#meaning if family then make sure the genera are different
										lowertaxa1 = species1.split(delims[i-1])[1].split(';')[0]
										lowertaxa2 = species2.split(delims[i-1])[1].split(';')[0]
										if lowertaxa1 != lowertaxa2:
											connected = 0
											for arg in argdict.keys():
												if ((species1 in argdict[arg][patient1]) and (species2 in argdict[arg][patient2])):
													connected = 1
											if connected == 1:
												linkedcount = linkedcount + 1
											allcount = allcount + 1
									else: #this is when it is genera and you already know that everything should be compared
										connected = 0
										for arg in argdict.keys():
											if ((species1 in argdict[arg][patient1]) and (species2 in argdict[arg][patient2])):
												connected = 1
										if connected == 1:
											linkedcount = linkedcount + 1
										allcount = allcount + 1
					if patient1 == patient2:
						withinbetween = 'within'
					else:
						withinbetween = 'between'
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(patient1, patient2, linkedcount, allcount, withinbetween,delimname))


