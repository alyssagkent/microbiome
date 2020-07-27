#gene distance to nearest linking hic-read
#so for every arg vs. taxonomy you will have a distance from that the start site of the ARG to the start site of the read (whatever is given in bwa output)
#to get there you have to track where the positions are of the reads, which comes from the hic trans reads table
#you could probably have a table that is scf1 scf2 a list of [(pos1, pos2), (pos1, pos2), pos1,pos2)], but you want to also include the reverse if you want to look it up
#so as you're iterating over the trans reads files, you put in for dict[scf1][scf2].append(pos1) but also dict[scf2][scf1].append(pos2)
#so then when you loop through the ARGs/taxa you have to keep track of extra information..
#another dictionary: [arg cluster][taxonomy] = [(ARG scf, taxa scf), (ARG scf, taxa scf)]
#because you're going to have multiple genes in each cluster, then you are going to find  
# then for each scaffold in that combination, check the distance between the ARG and the scaffold
#so you should probably make a function that does that cleanly
#so what are you going to need to plot this
#Sample Patient Cluster Taxonomy Minimum_Distance Distances
#I think it'd be useful to know how many scaffolds show this relationship?
#load prodigal information/dictionary from the simpler dictionary of their data

import glob
from Bio import SeqIO
from best_org import *
from igraph import *
import collections

depth = 1

#readdict given the argscf and the taxascf gives you a list of the positions of the reads linking the two on the argscf
readdict = {}
#check this path:
transpaths = glob.glob('/workdir/users/agk85/CDC/newhic/mapping/*/*_trans_primary_98.txt')
for transhandle in transpaths:
	with open(transhandle) as transfile:
		#header = transfile.readline() #does it have a header??
		for line in transfile:
			#check these positions
			scaffold1 = line.split('\t')[2]
			start1 = int(line.split('\t')[3])
			scaffold2 = line.split('\t')[8]
			start2 = int(line.split('\t')[9])
			try:
				readdict[scaffold1 + scaffold2].append(start1)
			except KeyError:
				readdict[scaffold1 + scaffold2] = [start1]
			try:
				readdict[scaffold2 + scaffold1].append(start2)
			except KeyError:
				readdict[scaffold2 + scaffold1]=[start2]


gene_pos_dict = {}
geneprothandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/arg_prot.fasta'
#change the handle if it is mobile genes to '/workdir/users/agk85/CDC/tables/metagenomes3/
with open(geneprothandle) as geneprotfile:
	for line in geneprotfile:
		if line[0] == '>':
			prot = line.split('>')[1].split(' ')[0]
			start = int(line.split(' # ')[1])
			stop = int(line.split(' # ')[2])
			gene_pos_dict[prot] = (start, stop)


def min_dist_reads(gene, taxascf):
	#this should get you the minimum distance between an ARG and a scaffold that has a taxonomy
	genescf = gene.split('_')[0] + '_' + gene.split('_')[1] + '_' + gene.split('_')[2]
	genepos_start, genepos_stop = gene_pos_dict[gene]
	readstarts = readdict[genescf + taxascf]
	readdists= []
	for readstart in readstarts:
		#get the minimum between the readstart and either the beginning or the end of the gene
		readstart_min = min(abs(readstart-genepos_start), abs(readstart-genepos_stop))
		#if the readstart is inside the gene
		if (genepos_start<readstart and readstart<genepos_stop):
			readstart_min = 0
		readdists.append(readstart_min)
	bestdist = min(readdists)
	return bestdist


def min_dist_scfs(gene, taxascfs):
	scaffolddists =[]
	for taxascf in taxascfs:
		scfdist = min_dist_reads(gene, taxascf)
		scaffolddists.append(scfdist)
	bestdist_fromscaffolds = min(scaffolddists)
	return bestdist_fromscaffolds



tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC' + '/combo_tables/metagenomes4/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()


best_org_dict = best_org(tbldict, depth)
g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping/trans_primary_ncol_' + '98' + '.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
#only keep the links greater than 2
g.es.select(weight_lt=2).delete()

samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

samples.sort()

inhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/args_95_nr.fna.clstr'
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
# argdict = collections.defaultdict(dict)
# for key in clusterprotdict.keys():
# 	argdict[key] = collections.defaultdict(dict)
# 	for samp in samples:
# 		argdict[key][samp] = [] #list for organisms

# argsampdict = {}
# for key in clusterprotdict.keys():
# 	argsampdict[key] = []
# 	genes = clusterprotdict[key]
# 	for gene in genes:
# 		gene_samp = gene.split('_')[0]
# 		argsampdict[key].append(gene_samp)

# samporgs = {}
# for samp in samples:
# 	samporgs[samp] = []

scflengths = {}
scfpaths = glob.glob("/workdir/users/agk85/CDC/prodigal_excise/metagenomes3/*/*_scaffold.fasta")
for scffile in scfpaths:
	for rec in SeqIO.parse(scffile, "fasta"):
		l = len(rec)
		scflengths[rec.id] = l

argdict = collections.defaultdict(dict)
#so basically as you are collecting the gene-organism connections you will want to populate a dictionary with the scaffolds that hold each taxonomy
#actually maybe you want to output each one when you loop through so you don't have to keep them---then you can min them by patient or by everyone
#so you will want to output this format
organisms = []
outdisthandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/min_dist/argdist_to_nearest_hicread.txt'
header = 'Cluster\tGene\tScaffold_length\tSample\tPatient\tMinimum_distance\tScf_taxonomy\tTaxonomy\n'
with open(outdisthandle,'w') as outfile:
	outfile.write(header)
	for cluster in clusterprotdict.keys():
		genes = clusterprotdict[cluster]
		for gene in genes:
			argdict[gene] = {}
			gene_samp = gene.split('_')[0]
			if gene_samp in samples:
				scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
				scaffold_length = scflengths[scf]
				try:
					arg_vids = g.neighborhood(scf, order=depth, mode=ALL)
					nodes_of_interest = g.vs[arg_vids]['name']
					for node in nodes_of_interest:
						besttaxa = best_org_dict[node]
						for besttaxon in besttaxa:
							if((besttaxon != '.') and ('Euk' not in besttaxon)):
								try:
									argdict[gene][besttaxon].append(node)
								except KeyError:
									argdict[gene][besttaxon] = [node]
					gene_taxonomies = argdict[gene].keys()
					if len(gene_taxonomies)>0:
						for taxa in gene_taxonomies:
							taxascfs = argdict[gene][taxa]
							if scf in taxascfs:
								in_genescf = '1'
								minimum_distance = 0
							else:
								in_genescf = '0'
								minimum_distance = min_dist_scfs(gene, taxascfs)
							patient = gene.split('-')[0]
							outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(cluster, gene, scaffold_length, gene_samp, patient, minimum_distance, in_genescf,taxa))
				except ValueError:
					a = 1 #nodes_of_interest = [scf]
