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
import collections
import sys


minreads = float(sys.argv[1])
genetype = sys.argv[2]

print('samples')
samples = []
samplist = glob.glob('/workdir/users/agk85/CDC2/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

samples.sort()



print('clusters')
if genetype == 'arg':
	inhandle = '/workdir/users/agk85/CDC2/args/args_99_nr.fna.clstr'
	geneprothandle = '/workdir/users/agk85/CDC2/args/arg_prot.fasta'

if genetype == 'mge':
	inhandle = '/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.fna.clstr'
	geneprothandle = '/workdir/users/agk85/CDC2/mobile/metagenomes/all_proteins.fna'


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



readdict = {}
contigdict = {}
scfdict = {}

print('hic_reads')

#check this path:
transpaths = glob.glob('/workdir/users/agk85/CDC2/hicpro/output/*/hic_results/data/*/*_allValidPairs')
for transhandle in transpaths:
	with open(transhandle) as transfile:
		#header = transfile.readline() #does it have a header??
		for line in transfile:
			#check these positions
			scaffold1 = line.split('\t')[1]
			start1 = int(line.split('\t')[2])
			scaffold2 = line.split('\t')[4]
			start2 = int(line.split('\t')[5])
			if scaffold1 != scaffold2:
				try:
					readdict[scaffold1 + scaffold2].append(start1)
				except KeyError:
					readdict[scaffold1 + scaffold2] = [start1]
				try:
					readdict[scaffold2 + scaffold1].append(start2)
				except KeyError:
					readdict[scaffold2 + scaffold1]=[start2]
				try:
					scfdict[scaffold1].append(scaffold2)
				except KeyError:
					scfdict[scaffold1]=[scaffold2]
				try:	
					scfdict[scaffold2].append(scaffold1)
				except KeyError:
					scfdict[scaffold2]=[scaffold1]


for scf in scfdict.keys():
	contigdict[scf] = list(set(scfdict[scf])) 

print('get gene positions')
gene_pos_dict = {}
#see above for geneprothandle
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


print('tbl dict')
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC2/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()

print('bin dict')
bindict = {}
bincontigs = glob.glob('/workdir/users/agk85/CDC2/das/*/*_DASTool_scaffolds2bin.txt')
for bincontig in bincontigs:
	with open(bincontig) as binfile:
		for line in binfile:
			binid = line.strip().split('\t')[1]
			contig = line.split('\t')[0]
			try:
				bindict[binid].append(contig)
			except KeyError:
				bindict[binid] = [contig]

print('contig vs bin')
genetaxa = {}
argdict = {}
arg_association={}
combocount = 0
combocount2 = 0
combocount3 = 0
combocount4=0
flagcount=0
cvbs = glob.glob('/workdir/users/agk85/CDC2/bins/*_das_1_contigs_v_bins.txt')
for cvb in cvbs:
	with open(cvb) as infile:
		header = infile.readline()
		for line in infile:
			quality = line.split('\t')[5]
			reads = float(line.split('\t')[8])
			if (quality != 'NA' and reads>=minreads):
				if genetype == 'arg':
					genes = line.split('\t')[16].split(',')
				if genetype == 'mge':
					genes = line.strip().split('\t')[19].split(',')
				taxonomy = line.split('\t')[12]
				if (taxonomy != 'k__; p__; c__; o__; f__; g__; s__;' and genes != ['NA'] and taxonomy != 'NA' and taxonomy != '.'):
					for gene in genes:
						combocount+=1
						contig1 = '_'.join(gene.split('_')[0:3])
						try:
							genetaxa[gene].append(taxonomy)
						except KeyError:
							genetaxa[gene] = [taxonomy]
						try:
							argdict[gene][taxonomy] = []
						except KeyError:
							argdict[gene] = {taxonomy:[]}
						#add all of the contigs that support 
						binid = line.split('\t')[3]
						combocount4+=1
						for contig2 in contigdict[contig1]:
							#if that contig is in the bin of interest
							if contig2 in bindict[binid]:
								argdict[gene][taxonomy].append(contig2)

print(combocount2, "fakers")

print('scflengths')
scflengths = {}
scfpaths = glob.glob("/workdir/users/agk85/CDC2/prodigal_excise/metagenomes/*/*_scaffold.fasta")
for scffile in scfpaths:
	for rec in SeqIO.parse(scffile, "fasta"):
		l = len(rec)
		scflengths[rec.id] = l


print('output and mindisting')
#so basically as you are collecting the gene-organism connections you will want to populate a dictionary with the scaffolds that hold each taxonomy
#actually maybe you want to output each one when you loop through so you don't have to keep them---then you can min them by patient or by everyone
#so you will want to output this format
organisms = []
counter = 0
failure_to_argdict= 0
emptytaxascfs = 0
outdisthandle = '/workdir/users/agk85/CDC2/args/min_dist/' + genetype + '_dist_to_nearest_hicread_min_' + str(minreads) + '.txt'
header = 'Cluster\tGene\tScaffold_length\tSample\tPatient\tMinimum_distance\tScf_taxonomy\tTaxonomy\tGenetype\tMinreads\n'
with open(outdisthandle,'w') as outfile:
	#output a header
	outfile.write(header)
	#iterate over the clusters
	for cluster in clusterprotdict.keys():
		genes = clusterprotdict[cluster]
		for gene in genes:
			scf = gene.split('_')[0] +'_' + gene.split('_')[1] +'_' + gene.split('_')[2]
			scaffold_length = scflengths[scf]
			try:
				gene_taxonomies = argdict[gene].keys()
				if len(gene_taxonomies)>0:
					for taxa in gene_taxonomies:
						taxascfs = argdict[gene][taxa]
						#this excludes phage repairing 0's essentially
						if (scf not in taxascfs and taxascfs != []):
							in_genescf = '0'
							minimum_distance = min_dist_scfs(gene, taxascfs)
							patient = gene.split('-')[0]
							gene_samp = '_'.join(gene.split('_')[0:3])
							outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(cluster, gene, scaffold_length, gene_samp, patient, minimum_distance, in_genescf,taxa,genetype, str(minreads)))
						else:
							emptytaxascfs +=1
			except KeyError:
				failure_to_argdict += 1


print("theoretical hic combos:", combocount)
print("actual unique hic combos:", combocount3)
print("failure to init argdict:", combocount2)
print("empty taxascfs:", emptytaxascfs)
print("failure to key argdict:", failure_to_argdict) 
