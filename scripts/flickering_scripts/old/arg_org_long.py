import numpy as np
import glob
import collections
from Bio import SeqIO
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-gc','--genecluster',dest='genecluster',action='store',required=True,type=str, help='Genes clstr file [Required]', metavar='ARGPID')
	parser.add_argument('-g','--genetype',dest='genetype',action='store',required=True,type=str, help='Gene [Required]', metavar='GENE')
	parser.add_argument('-m', '--minreads', dest='minreads', action='store', required=True, type=int, help='Min reads [Required]')
	parser.add_argument('-k', '--kraken_bins', dest='kraken_bins', action='store', required=True, help='Kraken_bins [Required]')
	parser.add_argument('-c','--checkm', dest="checkm", action='store', required=True, help='checkm_info [REQUIRED]', metavar="CHECKM")
	parser.add_argument('-i','--in', dest="cvb", action='store', required=True, help='Contig vs bin file [REQUIRED]', metavar="OUTFILE")
	parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="INFILE")
	args = parser.parse_args()
 	return(args)

args = getOptions()
####################################
#get args, get orgs cleanly, get arg-org connections based on a) taxa level, b) min reads (2 or 5)
#each will have a dictionary so you can check it, 
#genedict[sample]= list of genes of interest
#use what we have here
#orgdict[sample]= list of organisms passing anybin
#get from kraken and checkm
#geneorgdict[sample][gene][org] = 
#so get the args like normal
#so in the end you iterate over, every sample, every gene, every org, then ask if you connect

####BIN STUFF######
print('Checkm quality')
qualitydict = {}
cluster_length = {}
samplelist = []
with open(args.checkm) as checkm:
	for line in checkm:
		bin = line.split('\t')[0].split('.contigs')[0]
		sample = bin.split('.')[0].split('_')[0]
		samplelist.append(sample)
		completion = float(line.split('\t')[2])
		contamination = float(line.split('\t')[3])
		length = line.strip().split('\t')[5]
		quality = 'BAD'
		if (completion > 90 and contamination < 5):
			quality = 'HQ'
		elif (completion >= 50 and contamination < 10):
			quality = 'MQ'
		elif (completion < 50 and contamination < 10):
			quality = 'LQ'	
		qualitydict[bin] = quality
		cluster_length[bin] = length
		
samples = list(set(samplelist))
samples.sort()

print('Initialize orgfulldict')
orgfulldict = {}
orgdict = {}
for sample in samples:
	orgfulldict[sample] = []
	orgdict[sample] = []


organismsfull = []
with open(args.kraken_bins) as krakenfile:
	for line in krakenfile:
		bin = line.split('.contigs.fa.report.txt.besttaxid')[0]
		sample = bin.split('_')[0].split('.')[0]
		taxonomy = line.strip().split('\t')[7]
		orgfulldict[sample].append(taxonomy)
		organismsfull.append(taxonomy)

orgs = list(set(organismsfull))
orgs.sort()
for sample in samples:
	orgdict[sample] = list(set(orgfulldict[sample]))


print('Initialize the genefulldict')
genefulldict= {}
genedict = {}
for sample in samples:
	genefulldict[sample] = []
	genedict[sample] = []

with open(args.genecluster) as genefile:
	for line in genefile:
		cluster = line.split('\t')[0]
		genelist = line.split('\t')[2].split(',')
		for gene in genelist:
			sample = gene.split('_')[0]
			genefulldict[sample].append(cluster)

genesfull = []
for sample in samples:
	genedict[sample] = list(set(genefulldict[sample]))
	genesfull = genesfull + genedict[sample]

genes = list(set(genesfull))
genes.sort()

#count   sample  contig  bin     bin_length      quality association_type        total_count     trans_count     norm_l  norm_rf contig_taxonomy kraken_taxonomy gtdb_taxonomy   arg_presence    arg_clusters    arg_genes       mge_presence    mge_clusters    mge_genes
print('Initialize geneorgdict')
print(len(samples))
print(len(genes))
print(len(orgs))
good_quals = ['HQ','MQ','LQ']
geneorgdict = {}
for sample in samples:
	geneorgdict[sample] = {}
# 	for gene in genes:
# 		geneorgdict[sample][gene] = {}
# 		for org in orgs:
# 			geneorgdict[sample][gene][org]=0


print('Filling geneorgdict')
badorgs = []
badbins = []
allgenelist = []
with open(args.cvb) as hic:
	header = hic.readline()
	for line in hic:
		fields = line.strip().split('\t')
		trans_count = float(fields[8])
		quality = fields[5]
		association_type = fields[6]
		if args.genetype == 'arg':
			gene_pres = fields[14]
			genelist = fields[15].split(',')
		if args.genetype == 'mge':
			gene_pres = fields[17]
			genelist = fields[18].split(',')
		if ((quality in good_quals) and (gene_pres == '1')):
			if ((trans_count >= args.minreads) or (association_type=='cluster_resident')):
				sample = fields[1]
				contig = fields[2]
				binname = fields[3]
				taxonomy = fields[12]
				if (taxonomy != 'NA'):
					for gene in genelist:
						#print(gene, taxonomy)
						try:
							#print('gene1', gene)
							geneorgdict[sample][gene][taxonomy] = 1
							if taxonomy == 'k__; p__; c__; o__; f__; g__; s__;':
								geneorgdict[sample][gene][taxonomy] = 0
						except KeyError:
							#print('gene2', gene)
							allgenelist.append(gene)
							geneorgdict[sample][gene] = {}
							geneorgdict[sample][gene][taxonomy] = 1
							if taxonomy == 'k__; p__; c__; o__; f__; g__; s__;':
								geneorgdict[sample][gene][taxonomy] = 0

						

#only the genes of interest
genes = list(set(allgenelist))
print(len(genes))
#print(geneorgdict)


#print(geneorgdict['US8-1']['5381']['k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_sp._A1C1;'])

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
def get_code(gene, org, samp):
	genepres = '0'
	connection = '0'
	orgpres = '0'
	color = '6'
	try:
		if geneorgdict[samp][gene][org]==1:
			connection = '1'
	except KeyError:
		a = 1
	if gene in genedict[samp]:
		genepres = '1'
	if org in orgdict[samp]:
		orgpres = '1'
	if (genepres == '0') & (orgpres == '0') & (connection == '0'):
		color = '0'
	if (genepres == '0') & (orgpres == '1') & (connection == '0'):
		color = '1'
	if (genepres == '1') & (orgpres == '0') & (connection == '0'):
		color = '2'
	if (genepres == '1') & (orgpres == '1') & (connection == '0'):
		color = '3'
	if (genepres == '1') & (orgpres == '1') & (connection == '1'):
               color = '4'
	code = genepres + '\t' + orgpres + '\t' + connection + '\t' + color
	return code

header = 'Count\tCluster\tARG_name\tTop_ARG\tOrganism'
for samp in samples:
	header = header + '\t' + samp+'_gene\t' + samp + '_org\t' + samp + '_connect\t' + samp + 'color'

header = header +'\n'
count = 0
genecount = 0
toparg = 'NA'
argname = 'NA'
level = 'all'
print('Running through all the combinations')
outhandle = '/workdir/users/agk85/CDC2/bins/flickering/{0}_org_hic_separated_{1}_{2}.tbl'.format(args.genetype, level, str(args.minreads))
with open(outhandle,'w') as outfile:
	outfile.write(header)
	for gene in genes:
		#status
		genecount += 1
		if genecount % 50 == 0:
			print(float(genecount)/len(genes))
		#argname = arg_name_dict[clusternum_map[key]]
		for organism in orgs:
			count += 1
			tentativeline = '{0}\t{1}\t{2}\t{3}\t{4}'.format(str(count), gene, argname,toparg, organism)
			for samp in samples:
				code = get_code(gene, organism, samp)
				tentativeline = tentativeline +'\t' + code
			tentativeline = tentativeline + '\n'
			outfile.write(tentativeline)

