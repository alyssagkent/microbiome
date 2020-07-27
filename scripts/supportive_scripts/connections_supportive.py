import collections
from Bio import SeqIO
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organisms list given a genetype and minreads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-gc','--genecluster',dest='genecluster',action='store',required=True,type=str, help='Genes clstr file [Required]', metavar='ARGPID')
	parser.add_argument('-g','--genetype',dest='genetype',action='store',required=True,type=str, help='Gene [Required]', metavar='GENE')
	parser.add_argument('-m', '--minreads', dest='minreads', action='store', required=True, type=int, help='Min reads [Required]')
	parser.add_argument('-c','--checkm', dest="checkm", action='store', required=True, help='checkm_info [REQUIRED]', metavar="CHECKM")
	parser.add_argument('-i','--in', dest="cvb", action='store', required=True, help='Contig vs bin file [REQUIRED]', metavar="OUTFILE")
	parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="INFILE")
	args = parser.parse_args()
 	return(args)

args = getOptions()
####################################
#ok simpler you just want to know for each patient, out of all times when there is a connection, how consistent is it across the timepoints

####BIN STUFF######


samplelist = []
with open(args.binhandle) as binfile:
	for line in binfile:
		binid,patient,sample,quality,taxonomy = line.strip().split('\t')
		samplelist.append(sample)

samples = list(set(samplelist))
samples.sort()


print('Initialize the genefulldict')
genedict = {}
for sample in samples:
 	genedict[sample] = {}

with open(args.genecluster) as genefile:
	for line in genefile:
		cluster = line.split('\t')[0]
		genelist = line.split('\t')[2].split(',')
		for gene in genelist:
			sample = gene.split('_')[0]
			try:
				genedict[sample][cluster] = 1
			except KeyError:
				genedict[sample] = {}
				genedict[sample][cluster] = 1	

print('Initialize connection dict')
connection_dict = {}
patients = []
for sample in samples:
	connection_dict[sample] = {}

gt_assocs = {}
print('Filling connection dict')
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
		if ((quality !='NA') and (quality != 'BAD') and (gene_pres == '1')):
			if ((trans_count >= args.minreads) or (association_type=='cluster_resident_withouthicsupport') or (association_type=='cluster_resident_withhicsupport')):
				sample = fields[1]
				taxonomy = fields[12]
				if (taxonomy != 'NA' and taxonomy !='k__; p__; c__; o__; f__; g__; s__;'):
					for gene in genelist:
						#connection_dict[sample].append((gene,taxonomy,association_type))
						try:
							connection_dict[sample][gene+'%'+taxonomy].append(association_type)
						except KeyError:
							connection_dict[sample][gene+'%'+taxonomy] = [association_type]
#connect_dict = {}
#for sample in samples:
# 	connect_dict[sample] = list(set(connection_dict[sample]))


level = 'all'
outhandle = '/workdir/users/agk85/CDC2/bins_hicsupport/connections_{0}_org_{1}_{2}.txt'.format(args.genetype, level, str(args.minreads))
finished = {}
with open(outhandle,'w') as outfile:
	header = 'patient\tsample\tgenetype\tminthresh\tgeneid\tsource\ttaxonomy\n'
	outfile.write(header)
	for sample in samples:
		patient = sample.split('-')[0]
		gt_assocs = connection_dict[sample]
		for connection in gt_assocs.keys():
			geneid = connection.split('%')[0]
			taxonomy = connection.split('%')[1]
			associations = list(set(gt_assocs[connection]))
			if len(associations) ==1:
				source = associations[0]
			if len(associations)>1:	
				source = 'cluster_resident_withhicsupport'
			outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(patient,sample,args.genetype,args.minreads,geneid,source, taxonomy))

