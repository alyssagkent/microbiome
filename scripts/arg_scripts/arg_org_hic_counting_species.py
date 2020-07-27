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
hicpid = '99'
combohandle = '/workdir/users/agk85/' + folder + '/ComboDesign.txt'
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
	g = Graph.Read_Ncol('/workdir/users/agk85/' + folder + '/hicpro/single/' + hic + '_output/hic_results/data/' + hic + '/' + hic + '_trans_primary_' + hicpid + '_ncol_withexcise_noeuks.txt',weights=True, directed=True)
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
							#no core-core
							if (nodemge == '.' and scfmge == '.' and node!=scf):
								a = 1
							else:
								#no mobiletaxa to mobiletaxa
								if (nodemge == 'mge' and scfmge == 'mge' and node!=scf):
									a = 1
								else:
									for besttaxon in besttaxa:
										if besttaxon != '.':
											argdict[hic+mgm][key].append(besttaxon)
						except KeyError:
							besttaxa = ''
				except:
					a = 1

##############################################
print("onto org arg counting")
delims = ['s__','g__', 'f__','o__','c__','p__','k__']
delim_names=['species','genus','family','order','class','phylum','kingdom']
delim_names=['species']
#really you want the hic+mgm instead of pat
for i in range(len(delim_names)):
	header = 'Cluster\tARGname\tSample\tStudy\t' + delim_names[i] + '_count\tTaxa\n'
	arghisthandle = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/comparisons/nocorecorenomobilemobile/arg_orgcounts' + '_' + argpid + '_' + hicpid + '_' + delim_names[i] + str(depth) + '.txt'
	with open(arghisthandle, 'w') as argout:
		argout.write(header)
		delim = delims[i]
		print delim
		for combo in combos:
			hic = combo[0]
			mgm = combo[1]
			if 'B' in hic:
				study = 'CDC'
			if 'U' in hic:
				study = 'Healthy'
			arghistcount = {}
			for key in clusterprotdict.keys():
				patstring = ','.join(clusterpatientdict[key])
				if ('US' in patstring and 'B' in patstring):
					argname = arg_name_dict[clusternum_map[key]]
					arghistcount[key] = []
					for org in orglist:
						if org in argdict[hic+mgm][key]:
							level_org=org.split(delim)[1].split(';')[0]
							if level_org != '':
								arghistcount[key].append(level_org)
			for arg in clusterprotdict.keys():
				patstring = ','.join(clusterpatientdict[arg])
				if ('US' in patstring and 'B' in patstring):
					if hic.split('-')[0] in clusterpatientdict[arg]:
						argname = arg_name_dict[clusternum_map[arg]]
						s = list(set(arghistcount[arg]))
						l = len(s)		
						argout.write(arg + '\t' + argname + '\t' + hic + '\t' + study + '\t' + str(l) + '\t' + ','.join(s) +  '\n')
			#print hic
			#print arghistcount


