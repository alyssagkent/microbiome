#import numpy as np
#import matplotlib.pyplot as plt
import glob
import collections
import sys
from Bio import SeqIO
import igraph
#from igraph import *
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get mge versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
 	parser.add_argument('-d', '--depth', dest='depth', action='store', required=True, type=int, help='HiC depth [Required]', metavar='DEPTH')
	parser.add_argument('-f', '--folder', dest='folder', action='store', required=True, type=str, help='Folder [Required]', metavar='FOLDER')
	parser.add_argument('-c', '--contacts', dest='contact_thresh', action='store', required=True, type=str, help='Min contacts [Required]', metavar='THRESH') 
 	parser.add_argument('-b', '--cvb', dest='contig_v_bin', action='store', required=True, type=str, help='Contig_v_bin combined file [Required]', metavar='CONTACT_V_BIN')
	args = parser.parse_args()
 	return(args)

#make sure you run this command before
#cat ~/agk/CDC2/hicpro/single/*_output/hic_results/data/*/*_trans_primary_0_ncol_withexcise_noeuks.txt ~/agk/CDC2/hicpro/single/*_output/hic_results/data/*/*_trans_primary_0_ncol_withexcise_noeuks.txt > CDC+healthy_trans_primary_ncol_0_withexcise_noeuks.txt

 
args = getOptions()
folder = args.folder
depth = args.depth
mincontacts=args.contact_thresh

#tbldict
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/' + folder + '/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()

print(type(depth))

#so now the 0, 1, and 2 are 
best_org_dict = {}
with open(args.contig_v_bin) as infile:
	header = infile.readline()
	for line in infile:
		trans_count = float(line.split('\t')[8])
		if trans_count>=mincontacts:
			taxonomy=line.split('\t')[12]
			association_type=line.split('\t')[6]
			binofinterest=line.split('\t')[3]
			contig=line.split('\t')[2]
			if depth == 1 or depth == 2:
				if association_type == 'cluster_hic':
					try:
						best_org_dict[contig].append(taxonomy)
					except:
						best_org_dict[contig] = [taxonomy]
			if depth == 2:
				if association_type == 'contig_hic':
					try:
						best_org_dict[contig].append(taxonomy)
					except KeyError:
						best_org_dict[contig] = [taxonomy]


samples = []
samplist = glob.glob('/workdir/users/agk85/' + folder + '/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

#only keep the clusters that are machinery
machinerydict = {}
machinery_handle = '/workdir/users/agk85/' + folder + '/mobile/metagenomes/mge_99_nr.fna.clstr.desc'
with open(machinery_handle) as mach:
	header = mach.readline()
	for line in mach:
		cluster = line.split('\t')[0]
		phage_mach = line.split('\t')[3]
		machinerydict[cluster]=phage_mach

inhandle = '/workdir/users/agk85/' + folder + '/mobile/metagenomes/mge_' + '99' + '_nr.fna.clstr'
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

#this one uses 
for key in clusterprotdict.keys():
	genes = clusterprotdict[key]
	for gene in genes:
		contig = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
		try:	
			besttaxa = best_org_dict[contig]
			for besttaxon in besttaxa:
				try:
					if (besttaxon != '.' and besttaxon !='NA'):
						mgedict[key].append(besttaxon)
				except KeyError:
					continue
		except KeyError:
			continue

#initialize argdict

#mge machinery
descdict = {}
deschandle = '/workdir/users/agk85/' + folder + '/mobile/metagenomes/mge_99_nr.fna.clstr.desc'
with open(deschandle) as descfile:
	header = descfile.readline()
	for line in descfile:
		cluster = line.split('\t')[0]
		descriptions = line.split('\t')[1]
		descdict[cluster] = descriptions


sourcedict = {}
phagesourcedict = {}
sourcehandle = '/workdir/users/agk85/' + folder + '/mobile/metagenomes/mge_99_nr.fna.clstr.tbl.annot'
with open(sourcehandle) as sourcefile:
	header = sourcefile.readline()
	for line in sourcefile:
		cluster = line.split('\t')[0]
		source = line.split('\t')[3]
		phage = line.split('\t')[5]
		sourcedict[cluster] = source
		phagesourcedict[cluster] = phage
#######################separate org and arg and connection
protdict = SeqIO.to_dict(SeqIO.parse('/workdir/users/agk85/' + folder + '/mobile/metagenomes/mobile_proteins.fna', 'fasta'))

count = 0
protein_list = []
outhandle = '{0}/mobile/metagenomes/mge_machinery_base_hic_phage_withtaxa_{1}_{2}_{3}_{4}.tbl'.format('/workdir/users/agk85/' + folder , '99' , '99' ,str(args.depth), str(2))
with open(outhandle,'w') as outfile:
	outfile.write("Cluster\tGene\tSource\tNumber_base_taxa\tNumber_hic_taxa\tNumber_base_species\tNumber_hic_species\tBase_taxonomies\tHic_taxonomies\tBase_species\tHic_species\tPfams\n")
	for key in clusterprotdict.keys():
		pfams = descdict[key]
		source = sourcedict[key]
		phage = phagesourcedict[key]
		basetaxonomy_list='NA'
		hictaxonomy_list = mgedict[key]
		numbasetaxa = 0
		numhictaxa = len(set(hictaxonomy_list))
		if ((numbasetaxa >0 or numhictaxa>0) and int(phage)>0):
			hicspecies_list = []
			for taxonomy in hictaxonomy_list:
				taxa = taxonomy.split('s__')[1].split(';')[0]
				if taxa != '':
					hicspecies_list.append(taxa)
			hictaxonomies =  ','.join(list(set(hictaxonomy_list)))
			basetaxonomies = 'NA'
			hicspecies = ','.join(list(set(hicspecies_list)))
			basespecies = 'NA'
			numhicspecies = len(set(hicspecies_list))
			numbasespecies = 'NA'
			for gene in clusterprotdict[key]:
				rec = protdict[gene]
				protein_list.append(rec)
				outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(key, gene, source, numbasetaxa,numhictaxa, numbasespecies,numhicspecies, basetaxonomies,hictaxonomies, basespecies, hicspecies, pfams))

seqhandle = '{0}/mobile/metagenomes/mge_machinery_base_hic_phage_withtaxa_{1}_{2}_{3}_{4}.fna'.format('/workdir/users/agk85/'+folder , '99' , '99' , str(args.depth),str(2))
SeqIO.write(protein_list, seqhandle,"fasta")

