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
from best_org_mgm import *

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

g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping_mgm/trans_primary_ncol_' + '98' + '_withexcise_noeuks.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
g.es.select(weight_lt=2).delete()

#go through all the nodes, if excised phage, if the mother scaffold in the tbldict, add link weight = 2 between them in graph
#for node in tbldict.keys():
#	if 'phage' in node:
#		origscf = node.split('_')[0] + '_scaffold_' + node.split('_')[1].split('|')[0]
#		try:
#			a = tbldict[origscf]
#			try:
#				a = g.vs.find(node)
#			except ValueError:
#				#add vertex
#				g.add_vertex(node)
#			try:
#				a = g.vs.find(origscf)
#			except ValueError:
#				#add vertex
#				g.add_vertex(origscf)
#			g.add_edge(node, origscf, weight=2)
#		except KeyError:
#			continue

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
				arg_vids = g.neighborhood(scf, order=args.depth, mode=ALL)
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



#keycode 
# argpres | orgpres | connection
# 000 = 0
# 001 = Not possible
# 010 = 1
# 011 = Not possible
# 100 = 2
# 101 = Not possible
# 110 = 3
# 111 = 4

#######################separate org and arg and connection
def get_code(argdict, key, org, samp):
	argpres = '0'
	connection = '0'
	orgpres = '0'
	color = '6'
	if org in argdict[key][samp]:
		connection = '1'
	if samp in argsampdict[key]:
		argpres = '1'
	if org in sampleorgs[samp]:
		orgpres = '1'
	if (argpres == '0') & (orgpres == '0') & (connection == '0'):
		color = '0'
	if (argpres == '0') & (orgpres == '1') & (connection == '0'):
		color = '1'
	if (argpres == '1') & (orgpres == '0') & (connection == '0'):
		color = '2'
	if (argpres == '1') & (orgpres == '1') & (connection == '0'):
		color = '3'
	if (argpres == '1') & (orgpres == '1') & (connection == '1'):
               color = '4'
	code = argpres + '\t' + orgpres + '\t' + connection + '\t' + color
	return code

header = 'Count\tCluster\tARG_name\tTop_ARG\tOrganism'
for samp in samples:
	header = header + '\t' + samp+'_arg\t' + samp + '_org\t' + samp + '_connect\t' + samp + 'color'

header = header +'\n'

print("outputting combinations")

count = 0
#figure out if you want to keep the lines that only have the organism...about 1 million lines
outhandle = '{0}/arg_v_org/metagenomes3/mgm_histograms/arg_org_mgm_separated_{1}_{2}_{3}_{4}.tbl'.format('/workdir/users/agk85/CDC' , '95' , '98' ,str(args.depth), str(2))
with open(outhandle,'w') as outfile:
	outfile.write(header)
	for key in clusterprotdict.keys():
		argname = arg_name_dict[clusternum_map[key]]
		for organism in orglist:
			count += 1
			tentativeline = ''
			if key in top_args:
				toparg = cardres[key]
			else:
				toparg = 'NA'
			tentativeline = '{0}\t{1}\t{2}\t{3}\t{4}'.format(str(count), key, argname,toparg, organism)
			for samp in samples:
				code = get_code(argdict, key, organism, samp)
				tentativeline = tentativeline +'\t' + code
			tentativeline = tentativeline + '\n'
			outfile.write(tentativeline)




##############################################
##############################################
print("onto org arg counting")
delims = ['s__','g__', 'f__','o__','c__','p__','k__']
delim_names=['species','genus','family','order','class','phylum','kingdom']

patients = [samp.split('-')[0] for samp in samples]
patients = list(set(patients))
print(patients)
for pat in patients:
	for i in range(len(delims)):
		orghisthandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mgm_histograms/org_argcounts_'+ pat + '_' + delim_names[i] + str(args.depth) + '.txt'
		arghisthandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mgm_histograms/arg_orgcounts_'+ pat + '_' + delim_names[i] + str(args.depth) + '.txt'
		delim = delims[i]
		orghistcount = {}
		arghistcount = {}
		level_orgs = []
		for org in orglist:
			level_org = org.split(delim)[1].split(';')[0]
			if level_org !="":
				orghistcount[level_org] = []
				level_orgs.append(level_org)
		level_orgs_set = list(set(level_orgs))
		for key in clusterprotdict.keys():
			if pat in clusterpatientdict[key]:
				argname = arg_name_dict[clusternum_map[key]]
				arghistcount[key] = []
				for org in orglist:
					for samp in samples:
						patient = samp.split('-')[0]
						if patient == pat:
							if org in argdict[key][samp]:
								level_org=org.split(delim)[1].split(';')[0]
								if level_org != '':
									orghistcount[level_org].append(key)
									arghistcount[key].append(level_org)
		with open(orghisthandle,'w') as orgout:
			for level_org in level_orgs_set:
				s = list(set(orghistcount[level_org]))
				l = len(s)
				orgout.write(level_org + '\t' + str(l) + '\n')
		
		with open(arghisthandle, 'w') as argout:
			for arg in clusterprotdict.keys():
				if pat in clusterpatientdict[arg]:
					argname = arg_name_dict[clusternum_map[arg]]
					s = list(set(arghistcount[arg]))
					l = len(s)		
					argout.write(arg + '\t' + argname + '\t' + str(l) + '\n')


for i in range(len(delims)):
	orghisthandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mgm_histograms/org_argcounts_'+  delim_names[i] + str(args.depth) + '.txt'
	arghisthandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/mgm_histograms/arg_orgcounts_'+  delim_names[i] + str(args.depth) + '.txt'
	delim = delims[i]
	orghistcount = {}
	arghistcount = {}
	level_orgs = []
	for org in orglist:
		level_org = org.split(delim)[1].split(';')[0]
		if level_org !="":
			orghistcount[level_org] = []
			level_orgs.append(level_org)
	level_orgs_set = list(set(level_orgs))
	for key in clusterprotdict.keys():
		argname = arg_name_dict[clusternum_map[key]]
		arghistcount[key] = []
		for org in orglist:
			for samp in samples:
				if org in argdict[key][samp]:
					level_org=org.split(delim)[1].split(';')[0]						
					if level_org != '':									
						orghistcount[level_org].append(key)						
						arghistcount[key].append(level_org)						
	with open(orghisthandle,'w') as orgout:											    
		for level_org in level_orgs_set:											   
			s = list(set(orghistcount[level_org]))									     
			l = len(s)													 
			orgout.write(level_org + '\t' + str(l) + '\n')								     
																	   
	with open(arghisthandle, 'w') as argout:											   
		for arg in clusterprotdict.keys():											 
			argname = arg_name_dict[clusternum_map[arg]]								       
			s = list(set(arghistcount[arg]))										   
			l = len(s)													 
			argout.write(arg + '\t' + argname + '\t' + str(l) + '\n') 

