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
	description="""ARGs as they are associated with PATRIC database"""
	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-gc','--genecluster',dest='genecluster',action='store',required=True,type=str, help='Genes clstr file [Required]', metavar='ARGPID')
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

# inhandle = '/workdir/users/agk85/CDC2/args/args_99_nr.fna.clstr'
inhandle = args.genecluster
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

####read in gene-organism file and create a dictionary
organisms = []
# inhandle = '/workdir/users/agk85/CDC2/bins_hicsupport/connections_arg_org_all_2.txt' 
inhandle = args.connections
with open(inhandle) as infile:
	header = infile.readline()
	for line in infile:
		patient,sample,genetype,minreads,geneid,source,taxonomy = line.strip().split('\t')
		if source == 'cluster_resident_withhicsupport' or source =='cluster_resident_withouthicsupport':
			basedict[geneid].append(taxonomy)
		if source == 'cluster_hic':
			genedict[geneid].append(taxonomy)
		if (source == 'cluster_resident_withhicsupport' or source == 'cluster_hic'):
			genesupportdict[geneid].append(taxonomy)


namedict = {}
mechdict = {}
broadmechdict = {}
deschandle = '/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech.txt'
with open(deschandle) as descfile:
	header = descfile.readline()
	for line in descfile:
		cluster = line.split('\t')[0]
		descriptions = line.split('\t')[2]
		mechanism = line.split('\t')[5]
		sub_mechanism=line.strip().split('\t')[6]
		namedict[cluster] = descriptions
		mechdict[cluster] = mechanism
		broadmechdict[cluster] = sub_mechanism


#######################separate org and arg and connection
def get_code(genedict, key, org):
	connection = '0'
	if org in genedict[key]:
		connection = '1'
	return connection

print("outputting combinations")

#protdict = SeqIO.to_dict(SeqIO.parse("/workdir/users/agk85/CDC/arg_v_org/metagenomes3/arg_prot.fasta", "fasta"))
protdict = SeqIO.to_dict(SeqIO.parse("/workdir/users/agk85/CDC2/args/args_99_nr.fna","fasta"))
count = 0
protein_list = []
outhandle = args.outhandle
with open(outhandle,'w') as outfile:
	outfile.write("Cluster\tGene\tNumber_genes_in_cluster\tName\tBroad_mechanism\tNumber_base_taxa\tNumber_hic_taxa\tBase_taxonomies\tHic_taxonomies\tHic_taxonomies_supportive\n")
	for key in clusterprotdict.keys():
		name = namedict[key]
		broadmechanism = broadmechdict[key]
		basetaxonomy_list=basedict[key]
		hictaxonomy_list = genedict[key]
		hic_supportive_taxonomy_list = genesupportdict[key]
		numbasetaxa = len(set(basetaxonomy_list))
		numhictaxa = len(set(hictaxonomy_list))
		numhicsupporttaxa = len(set(hic_supportive_taxonomy_list))
		if (numbasetaxa >0 or numhictaxa>0 or numhicsupporttaxa>0):
			hictaxonomies =  ','.join(list(set(hictaxonomy_list)))
			basetaxonomies = ','.join(list(set(basetaxonomy_list)))
			hicsupportivetaxonomies =  ','.join(list(set(hic_supportive_taxonomy_list)))
			gene = clusternum_map[key]
			numgenes = str(len(clusterprotdict[key]))
			rec = protdict[gene]
			protein_list.append(rec)
			outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(key, gene,numgenes, name,broadmechanism, numbasetaxa,numhictaxa, basetaxonomies,hictaxonomies,numhicsupporttaxa, hicsupportivetaxonomies))


# seqhandle = '{0}/args/patric_comparisons/arg_base_hic_withtaxa_{1}_{2}_{3}_{4}.fna'.format('/workdir/users/agk85/CDC2' , '99' , '99' ,str(2))
# seqhandle = args.outhandle
# SeqIO.write(protein_list, seqhandle,"fasta")

