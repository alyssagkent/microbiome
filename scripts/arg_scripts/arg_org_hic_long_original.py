#refined ARG vs. Organism traversal
# #use the best organism function from CDC one
#python arg_org_hic_refined.py -b /workdir/users/agk/gates -a 95 -p 0 -d 2
#goal of latest update
#instead of using just the best organism given by the combo_tables, instead use the best organism given by the clusters from a given cluster file
import numpy as np
#import matplotlib.pyplot as plt
import glob
import collections
#import sys
from Bio import SeqIO
from igraph import *
import argparse
from argparse import RawDescriptionHelpFormatter
from best_org_folder import *

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-a','--argpid',dest='argpid',action='store',required=True,type=str, help='ARG clusterd at PID [Required]', metavar='ARGPID')
	parser.add_argument('-p','--hicpid',dest='hicpid',action='store',required=True,type=str, help='HIC mapped at PID [Required]', metavar='HICPID')
	parser.add_argument('-f', '--folder', dest='folder', action='store', required=True, type=str, help='Folder [Required]', metavar='FOLDER')
 	parser.add_argument('-d', '--depth', dest='depth', action='store', required=True, type=int, help='HiC depth [Required]', metavar='DEPTH') 
 	parser.add_argument('-c', '--combinations', dest='combofile', action='store', required=True, type=str, help='Combination file [Required]', metavar='COMBINATIONS')
	args = parser.parse_args()
 	return(args)

#for testing
depth = 1
folder = 'press'
argpid = '99'
hicpid = '98'
combohandle = '/workdir/users/agk85/press/ComboDesign.txt'
minreads=2
cb1 = 'MluCI-1ProxiMeta-1'
cb2 = 'Sau3aI-1ProxiMeta-1'


args = getOptions()
depth = args.depth
folder = args.folder
argpid = args.argpid
hicpid = args.hicpid
combohandle = args.combofile
minreads = 2

#tbldict
tbldict = {}
mgedict = {}
tblpaths = glob.glob('/workdir/users/agk85/' + folder + '/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()
			mgedict[line.split('\t')[0]] = line.strip().split('\t')[-2]


#best_org_dict[hic+mgm] = {}
#best_org_d = best_org(tbldict, depth, folder, hicpid, minreads)

#import the hic to metagenome connections as tuples
combos = []
with open(combohandle) as cbfile:
	for line in cbfile:
		hic = line.split('\t')[0]
		mgm = line.split('\t')[1].strip()
		combos.append((hic, mgm))


best_org_dict = {}
for combo in combos:
	hic = combo[0]
	mgm = combo[1]
	cb = hic+mgm
	print hic, mgm
	best_org_d = best_org_folder(tbldict, depth, folder, hicpid, minreads, hic, mgm)	
	best_org_dict[cb] = best_org_d

#legacy keep this to maintain the proper spacing
top_args = []

#make a dictionary of graphs based on the hic library it uses
#get the library from the name somehow
#g = Graph.Read_Ncol('/workdir/users/agk85/' + folder + '/hic/mapping/trans_primary_ncol_' + hicpid + '_withexcise_noeuks.txt',weights=True, directed=True)
#g.to_undirected(mode="collapse",combine_edges="sum")
#g.es.select(weight_lt=minreads).delete()
	
#graphs
graphs = {}
for combo in combos:
	hic = combo[0]
	mgm = combo[1]
	g = Graph.Read_Ncol('/workdir/users/agk85/' + folder + '/hic/mapping/' + hic + '/' + hic + '_' + mgm+ '_trans_primary_' + hicpid + '_ncol_withexcise_noeuks.txt',weights=True, directed=True)
	g.to_undirected(mode="collapse",combine_edges="sum")
	g.es.select(weight_lt=minreads).delete()
	graphs[hic+mgm] = g 

samples = []
samplist = glob.glob('/workdir/users/agk85/' + folder + '/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

inhandle = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/args_' + argpid + '_nr.fna.clstr'
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
prot_arg_name = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/mapping/bwa_alignments_' + argpid +'_'+ argpid+'/arg_v_samp_'+argpid+ '_' + argpid + '_names.txt'
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
				#this means that the protein is not in the set of samples we are interested in 

#make a shortened list of argnames---i don't think i use this ever.
#argnames = list(set(arg_name_dict.values()))
#argnames.sort()


#initialize argdict
argdict = collections.defaultdict(dict)
argsampdict = {}
for combo in combos:
	cb = combo[0] + combo[1]
	argdict[cb] = collections.defaultdict(dict)
	argsampdict[cb] = collections.defaultdict(dict)
	for key in clusterprotdict.keys():
		argdict[cb][key] = [] #list for organisms
		argsampdict[cb][key] = []
		genes = clusterprotdict[key]
		for gene in genes:
			gene_samp = gene.split('_')[0]
			argsampdict[cb][key].append(gene_samp)


#do i need a cb component?
samporgs = {}
sampleorgs = {}
organisms = []
for combo in combos:
	cb = combo[0] + combo[1]
	samporgs[cb] = []
	for node in best_org_dict[cb]:
		besttaxa = best_org_dict[cb][node]
		for besttaxon in besttaxa:
			if besttaxon != '.':
				samporgs[cb].append(besttaxon)
				organisms.append(besttaxon)
        sampleorgs[cb] = list(set(samporgs[cb]))

orglist = list(set(organisms))
orglist.sort()

#modified to separate hic and mgm libraries
for combo in combos:
	hic = combo[0]
	mgm = combo[1]
	cb = hic+mgm
	g = graphs[cb]
	for key in clusterprotdict.keys():
		genes = clusterprotdict[key]
		for gene in genes:
			gene_samp = gene.split('_')[0]
			if gene_samp==mgm:
				scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
				try:
					scfmge = mgedict[scf] #this determines the scaffolds mobility
					try:
						arg_vids = g.neighborhood(scf, order=depth, mode=ALL)
						nodes_of_interest = g.vs[arg_vids]['name']
					except ValueError:
						#this is when the scaffold isn't in the network
						nodes_of_interest = [scf]
					for node in nodes_of_interest:
						try:
							besttaxa = best_org_dict[cb][node]
							nodemge = mgedict[node]
							for besttaxon in besttaxa:
								if besttaxon != '.':
									argdict[hic+mgm][key].append(besttaxon)
						except KeyError:
							besttaxa = ''
				except:
					a = 1

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
def get_code(argdict,hic, mgm, key, org):
	argpres = '0'
	connection = '0'
	orgpres = '0'
	color = '6'
	cb = hic+mgm
	if org in argdict[cb][key]:
		connection = '1'
	if mgm in argsampdict[cb][key]:
		argpres = '1'
	if org in sampleorgs[cb]:
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
for combo in combos:
	hic = combo[0]
	mgm = combo[1]
	header = header + '\t' + hic+ '+' + mgm +'_arg\t' + hic+ '+' + mgm  + '_org\t' + hic+ '+' + mgm  + '_connect\t' + hic+ '+' + mgm + 'color'

header = header +'\n'

print("outputting combinations")

count = 0
#figure out if you want to keep the lines that only have the organism...about 1 million lines
outhandle = '{0}/arg_v_org/metagenomes/histograms/original/arg_org_hic_separated_{1}_{2}_{3}_{4}.tbl'.format('/workdir/users/agk85/' + folder , argpid , hicpid ,str(depth), str(2))
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
			for combo in combos:
				hic = combo[0]
				mgm = combo[1]
				code = get_code(argdict, hic, mgm, key, organism)
				tentativeline = tentativeline +'\t' + code
			tentativeline = tentativeline + '\n'
			outfile.write(tentativeline)
##############################################
##############################################
print("onto org arg counting")
delims = ['s__','g__', 'f__','o__','c__','p__','k__']
delim_names=['species','genus','family','order','class','phylum','kingdom']

#really you want the hic+mgm instead of pat
for combo in combos:
	hic = combo[0]
	mgm = combo[1]
	for i in range(len(delims)):
		arghisthandle = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/histograms/original/arg_orgcounts_'+ hic + '+' + mgm + '_' + argpid + '_' + hicpid + '_' + delim_names[i] + str(depth) + '.txt'
		delim = delims[i]
		arghistcount = {}
		level_orgs = []
		for org in orglist:
			level_org = org.split(delim)[1].split(';')[0]
			if level_org !="":
				level_orgs.append(level_org)
		level_orgs_set = list(set(level_orgs))
		for key in clusterprotdict.keys():
			argname = arg_name_dict[clusternum_map[key]]
			arghistcount[key] = []
			for org in orglist:
				if org in argdict[hic+mgm][key]:
					level_org=org.split(delim)[1].split(';')[0]
					if level_org != '':
						arghistcount[key].append(level_org)
		with open(arghisthandle, 'w') as argout:
			for arg in clusterprotdict.keys():
				argname = arg_name_dict[clusternum_map[arg]]
				s = list(set(arghistcount[arg]))
				l = len(s)		
				argout.write(arg + '\t' + argname + '\t' + str(l) + '\n')

for i in range(len(delims)):
	arghisthandle = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/histograms/original/arg_orgcounts_'+ argpid + '_' + hicpid + '_' + delim_names[i] + str(depth) + '.txt'
	delim = delims[i]
	arghistcount = {}
	level_orgs = []
	for org in orglist:
		level_org = org.split(delim)[1].split(';')[0]
		if level_org !="":
			level_orgs.append(level_org)
	level_orgs_set = list(set(level_orgs))
	for key in clusterprotdict.keys():
		argname = arg_name_dict[clusternum_map[key]]
		arghistcount[key] = []
		for org in orglist:
			for combo in combos:
				cb = combo[0] + combo[1]
				if org in argdict[cb][key]:
					level_org=org.split(delim)[1].split(';')[0]
					if level_org != '':
						arghistcount[key].append(level_org)
	with open(arghisthandle, 'w') as argout:
		for arg in clusterprotdict.keys():
			argname = arg_name_dict[clusternum_map[arg]]
			s = list(set(arghistcount[arg]))
			l = len(s)
			argout.write(arg + '\t' + argname + '\t' + str(l) + '\n') 

