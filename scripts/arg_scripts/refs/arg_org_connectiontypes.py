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
from best_org import best_org

#combo tables
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC' + '/combo_tables/metagenomes4/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()

depth = 1
best_org_dict = best_org(tbldict, depth)

g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping/trans_primary_ncol_' + '98' + '_noeuks.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
g.es.select(weight_lt=2).delete()


typedict = {}
for scf in tbldict:
	flag = 0
	mobility = tbldict[scf].split('\t')[-2]
	if mobility == 'mge':
		typedict[scf] = 'mobile'
		flag = 1
	if mobility == '.':
		besttaxa = best_org_dict[scf]
		taxa = []
		for besttaxon in besttaxa:
			if besttaxon != '.':
				taxa.append(besttaxon)
		if len(taxa)==0:
			typedict[scf] = 'nonmobilenotaxa'
			flag = 1
		if len(taxa) >0:
			typedict[scf] = 'genomic'
			flag =1
	if flag == 0:
		print scf

samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
        samp = samplong.split('/')[6]
        samples.append(samp)

patients = ['B314','B316','B320','B331','B335','B357','B370']
gs = ['g1','g2','g3','g4','g5','g6']
#g1-mobile mobile, g2-mobile nonmobilenotaxa,g3-mobile genomic,g4-nonmobilenotaxa nonmobilenotaxa,g5-nonmobilenotaxa genomic,g6-genomic genomic
groups = {}
for pat in patients:
	groups[pat] = {}
	for gnum in gs:
		groups[pat][gnum] = 0
#groups 
for edge in g.es:
	scf1 = g.vs[edge.source]['name']
	scf2 = g.vs[edge.target]['name']
	patient = scf1.split('-')[0]
	type1 = typedict[scf1]
	type2 = typedict[scf2]
	if type1 == 'mobile':
		if type2 == 'mobile':
			groups[patient]['g1'] += 1
		if type2 == 'nonmobilenotaxa':
			groups[patient]['g2']  += 1
		if type2 == 'genomic':
			groups[patient]['g3']  += 1
	if type1 == 'nonmobilenotaxa':
		if type2 == 'mobile':
			groups[patient]['g2']  += 1
		if type2 == 'nonmobilenotaxa':
			groups[patient]['g4']  += 1
		if type2 == 'genomic':
			groups[patient]['g5']  += 1	
	if type1 == 'genomic':	
		if type2 == 'mobile':
			groups[patient]['g3']  += 1
		if type2 == 'nonmobilenotaxa':
			groups[patient]['g5']  += 1
		if type2 == 'genomic':
			groups[patient]['g6']  += 1

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


genegroups = {}
for pat in patients:
        genegroups[pat] = {}
        for gnum in gs:
                genegroups[pat][gnum] = 0

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
								patient = scf.split('-')[0]
								type1 = typedict[scf]
								type2 = typedict[node]
								if type1 == 'mobile':
									if type2 == 'mobile':
										genegroups[patient]['g1'] += readcount
									if type2 == 'nonmobilenotaxa':
										genegroups[patient]['g2']  += readcount
									if type2 == 'genomic':
										genegroups[patient]['g3']  += readcount
								if type1 == 'nonmobilenotaxa':
									if type2 == 'mobile':
										genegroups[patient]['g2']  += readcount
									if type2 == 'nonmobilenotaxa':
										genegroups[patient]['g4']  += readcount
									if type2 == 'genomic':
										genegroups[patient]['g5']  += readcount
								if type1 == 'genomic': 
									if type2 == 'mobile':
										genegroups[patient]['g3']  += readcount
									if type2 == 'nonmobilenotaxa':
										genegroups[patient]['g5']  += readcount
									if type2 == 'genomic':
										genegroups[patient]['g6']  += readcount
								try:
									argdict[key][gene_samp][besttaxon]+= readcount
								except KeyError:
									argdict[key][gene_samp][besttaxon] = readcount
					except KeyError:
						besttaxa = ''
			except ValueError:
				a=1


outhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/rarefaction/groups.txt'
with open(outhandle,'w') as outfile:
	outfile.write("Patient\t" + '\t'.join(gs)+'\n')
	for patient in patients:
		outfile.write(patient)
		for g in gs:
			outfile.write('\t' + str(groups[patient][g]))
		outfile.write('\n')
		

outhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/rarefaction/groups_linkingARGs.txt'
with open(outhandle,'w') as outfile:
        outfile.write("Patient\t" + '\t'.join(gs)+'\n')
	for patient in patients:
                outfile.write(patient)
                for g in gs:
                        outfile.write('\t' + str(genegroups[patient][g]))
                outfile.write('\n')
