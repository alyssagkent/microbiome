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
orgdict = {}
with open(args.binhandle) as binfile:
	for line in binfile:
		binid,patient,sample,quality,taxonomy = line.strip().split('\t')
		samplelist.append(sample)
		try:
			a = orgdict[sample]
		except KeyError:
			orgdict[sample] = {}
		if taxonomy != '.':
			orgdict[sample][taxonomy] = 1

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



print('Initialize geneorgdict')
geneorgdict = {}
connection_dict = {}
patients = []
for sample in samples:
	geneorgdict[sample] = {}
	connection_dict[sample] = []
	patient = sample.split('-')[0]
	patients.append(patient)


patients = list(set(patients))
patsamps = {}
for patient in patients:
	patsamps[patient] = []

for sample in samples:
	pat = sample.split('-')[0]
	patsamps[pat].append(sample)

print('Filling geneorgdict')
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
		if ((quality !='BAD') and (gene_pres == '1')):
			if ((trans_count >= args.minreads) or (association_type=='cluster_resident')):
				sample = fields[1]
				taxonomy = fields[12]
				if (taxonomy != 'NA' and taxonomy !='k__; p__; c__; o__; f__; g__; s__;'):
					for gene in genelist:
						try:
							geneorgdict[sample][gene][taxonomy] = 1
						except KeyError:
							geneorgdict[sample][gene] = {taxonomy:1}
						connection_dict[sample].append((gene,taxonomy))
						

level = 'all'
print('Running through all the combinations')
outhandle = '/workdir/users/agk85/CDC2/bins/flickering2/{0}_org_flickering_{1}_{2}.txt'.format(args.genetype, level, str(args.minreads))
connected_dict = {}
unconnected_dict = {}
patient_connection_dict = {}
for patient in patients:
	connected_dict[patient] = 0
	unconnected_dict[patient] = 0
	patient_connection_dict[patient] = []

patient_connection_dict = {}
for sample in samples:
	print(len(connection_dict[sample]))
	patient = sample.split('-')[0]
	

for patient in patients:
	print(len(connection_dict[sample]))



assessed = []
for sample in samples:
	print(sample)
	connections = connection_dict[sample]
	#get the patient	
	patient = sample.split('-')[0]
	patientsamples = patsamps[patient]
	for connection in connections:
		gene = connection[0]
		taxonomy = connection[1]
		#iterate through the samples
		#is the gene present
		#is the taxa present?
		#are they connected?
		#if 1 and 2 but not 3 update patient unconnected total
		#if 1 and 2 and 3 update patient connected total
		#you don't want to overcount 
		if (sample,gene,taxonomy) not in assessed:
			for samp in patientsamples:
				orgpres = 0
				genepres = 0
				connected = 0
				#print(samp,gene,taxonomy)
				try:
					orgpres = orgdict[samp][taxonomy]
					#print(samp,taxonomy)
				except KeyError:
					orgpres = 0
					#print('no org')
				try:
					genepres = genedict[samp][gene]
					#print(samp,gene)
				except KeyError:
					genepres = 0
					#print('no gene')
				try: 
					connected = geneorgdict[samp][gene][taxonomy]
					assessed.append((samp,gene,taxonomy))
					#print('connected')
				except KeyError:
					connected = 0
					#print('unconnected')
				if (orgpres and genepres and not connected):
					unconnected_dict[patient] += 1
					#print('report_un')
				if orgpres and genepres and connected:
					connected_dict[patient] += 1
					#print('report_con')

with open(outhandle,'w') as outfile:
	header = 'patient\tgenetype\tminthresh\tconsistency\n'
	outfile.write(header)
	for patient in patients:
		print(unconnected_dict[patient])
		print(connected_dict[patient])
		consistency = str(float(connected_dict[patient])/float(connected_dict[patient]+unconnected_dict[patient]))
		outfile.write('{0}\t{1}\t{2}\t{3}\n'.format(patient,args.genetype,args.minreads,consistency))

