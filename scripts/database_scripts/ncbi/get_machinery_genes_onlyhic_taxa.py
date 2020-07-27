#refined ARG vs. Organism traversal
# #use the best organism function from CDC one
#python arg_org_hic_refined.py -b /workdir/users/agk/gates -a 95 -p 0 -d 2
#goal of latest update
#instead of using just the best organism given by the combo_tables, instead use the best organism given by the clusters from a given cluster file
import glob
import collections
import sys
from Bio import SeqIO
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-gc','--genecluster',dest='genecluster',action='store',required=True,type=str, help='Genes clstr file [Required]', metavar='ARGPID')
	parser.add_argument('-m',dest='minreads',action='store',required=True,type=int, help='Minreads [Required]', metavar='MINREADS')
	parser.add_argument('-o',dest='outhandle',action='store',required=True,type=str, help='Outfile [Required]', metavar='OUTFILE')
	parser.add_argument('-c','--connections', dest="connections", action='store', required=True, help='Contig vs bin file [REQUIRED]', metavar="CONNECTIONS")
	args = parser.parse_args()
 	return(args)
 
args = getOptions()

#tbldict
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC2/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()


samples = []
samplist = glob.glob('/workdir/users/agk85/CDC2/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[-1]
	samples.append(samp)

#only keep the clusters that are machinery
machinerydict = {}
machinery_handle = '/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.fna.clstr.desc'
with open(machinery_handle) as mach:
	header = mach.readline()
	for line in mach:
		cluster = line.split('\t')[0]
		phage_mach = line.split('\t')[3]
		machinerydict[cluster]=phage_mach

inhandle = '/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.fna.clstr'
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


#initialize genedict
genedict = collections.defaultdict(dict)
genesupportdict = collections.defaultdict(dict)
basedict = collections.defaultdict(dict)
for key in clusterprotdict.keys():
	genedict[key] = []
	basedict[key] = []
	genesupportdict[key] = []


# inhandle = '/workdir/users/agk85/CDC2/bins_hicsupport/connections_arg_org_all_2.txt'
inhandle = args.connections
with open(inhandle) as infile:
	header = infile.readline()
	for line in infile:
		patient,sample,genetype,minreads,geneid,source,taxonomy = line.strip().split('\t')
		try:
			empty = clusterprotdict[geneid]
			if source == 'cluster_resident_withhicsupport' or source =='cluster_resident_withouthicsupport':
				basedict[geneid].append(taxonomy)
			if source == 'cluster_hic':
				genedict[geneid].append(taxonomy)
			if (source == 'cluster_resident_withhicsupport' or source=='cluster_hic'):
				genesupportdict[geneid].append(taxonomy)
		except KeyError:
			#the gene is not a phage machinery
			a = 1
	
#mge machinery
descdict = {}
deschandle = '/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.fna.clstr.desc'
with open(deschandle) as descfile:
	header = descfile.readline()
	for line in descfile:
		cluster = line.split('\t')[0]
		descriptions = line.split('\t')[1]
		descdict[cluster] = descriptions


sourcedict = {}
phagesourcedict = {}
sourcehandle = '/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.fna.clstr.tbl.annot'
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

protdict = SeqIO.to_dict(SeqIO.parse("/workdir/users/agk85/CDC2/mobile/metagenomes/mobile_proteins.fna", "fasta"))

count = 0
protein_list = []
outhandle = args.outhandle
with open(outhandle,'w') as outfile:
	#from arg replaced name with source
	#from arg replace Broad_mechanism with pfams
	outfile.write("Cluster\tGene\tNumber_genes_in_cluster\tSource\tPfams\tNumber_base_taxa\tNumber_hic_taxa\tBase_taxonomies\tHic_taxonomies\tNumber_hicsupport_taxa\tHic_taxonomies_supportive\n")
	#outfile.write("Cluster\tGene\tSource\tNumber_base_taxa\tNumber_hic_taxa\tNumber_base_species\tNumber_hic_species\tBase_taxonomies\tHic_taxonomies\tBase_species\tHic_species\tPfams\n")
	for key in clusterprotdict.keys():
		pfams = descdict[key]
		source = sourcedict[key]
		phage = phagesourcedict[key]
		basetaxonomy_list=basedict[key]
		hictaxonomy_list = genedict[key]
		hicsupporttaxonomy_list = genesupportdict[key]	
		numbasetaxa = len(set(basetaxonomy_list))
		numhictaxa = len(set(hictaxonomy_list))
		numhicsupporttaxa = len(set(hicsupporttaxonomy_list))
		if ((numbasetaxa >0 or numhictaxa>0 or numhicsupporttaxa>0) and int(phage)>0):
			hictaxonomies =  ','.join(list(set(hictaxonomy_list)))
			basetaxonomies = ','.join(list(set(basetaxonomy_list)))
			hicsupporttaxonomies=','.join(list(set(hicsupporttaxonomy_list)))
			numgenes = str(len(clusterprotdict[key]))
			gene = clusternum_map[key]
			outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(key, gene, numgenes, source, pfams, numbasetaxa,numhictaxa,basetaxonomies,hictaxonomies,numhicsupporttaxa,hicsupporttaxonomies))

seqhandle = '{0}/mobile/metagenomes/mge_machinery_base_hic_phage_withtaxa_{1}_{2}.fna'.format('/workdir/users/agk85/CDC2' , '99' , str(args.minreads))
SeqIO.write(protein_list, seqhandle,"fasta")

