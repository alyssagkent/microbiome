#smillie esque taxa comparisons at genus level first
#improve to looking at different taxonomic levels
import numpy as np
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
patients = []
samplist = glob.glob('/workdir/users/agk85/CDC2/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)
	patients.append(samp.split('-')[0])

patients = list(set(patients))
patients.sort()

inhandle = args.genecluster
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
argdict = collections.defaultdict(dict)
for key in clusterprotdict.keys():
	argdict[key] = collections.defaultdict(dict)
	for pat in patients:
		argdict[key][pat] = [] #list for organisms

#initalize argpatdict
argpatdict = {}
for key in clusterprotdict.keys():
	argpatdict[key] = []
	genes = clusterprotdict[key]
	for gene in genes:
		gene_pat = gene.split('-')[0]
		argpatdict[key].append(gene_pat)

#initialize patorgs
patorgs = {}
for pat in patients:
	patorgs[pat] = []


####read in gene-organism file and create a dictionary
organisms = []
with open(args.connections) as infile:
	header = infile.readline()
	for line in infile:
		patient,sample,genetype,minreads,geneid,taxonomy = line.strip().split('\t')
		argdict[geneid][patient].append(taxonomy) 
		organisms.append(taxonomy)
		patorgs[patient].append(taxonomy)


#make lists unique
patientorgs = {}
for pat in patorgs:
	patientorgs[pat] = list(set(patorgs[pat]))


#get all of the organisms...make a set
orglist = list(set(organisms))
orglist.sort()

#################################################################ok so let's begin
#taxalists
delims = ['g__','f__','o__','c__','p__','k__']
delimnames = ['genus','family','order','class','phylum','kingdom']
taxalists = collections.defaultdict(dict)
for patient in patients:
	taxalists[patient] = collections.defaultdict(dict)
	for delimname in delimnames:
		taxalists[patient][delimname]= collections.defaultdict(dict)


alltaxa = collections.defaultdict(dict)

for i in range(len(delimnames)):
	delimname = delimnames[i]
	delim = delims[i]
	alltaxa[delimname] = []
	for patient in patients:
		for org in patientorgs[patient]:
			if org.split('s__')[1].split(';')[0] != '':
				taxa = org.split(delim)[1].split(';')[0]
				if taxa != '':
					alltaxa[delimname].append(taxa)
					try:
						taxalists[patient][delimname][taxa].append(org)
					except AttributeError:
						taxalists[patient][delimname][taxa] = [org]


for delimname in delimnames:
	alltaxa[delimname] = list(set(alltaxa[delimname]))
	alltaxa[delimname].sort()

#comparing two people
finished = []
outhandle = args.outhandle 
with open(outhandle, 'w') as outfile:
	for patient1 in patients:
		for patient2 in patients:
			if (patient2,patient1) not in finished:
				finished.append((patient1,patient2))
				for i in range(len(delimnames)):
					delimname = delimnames[i]
					linkedcount = 0
					allcount = 0
					for taxa in alltaxa[delimname]:
						t1 = taxa in taxalists[patient1][delimname].keys()
						t2 = taxa in taxalists[patient2][delimname].keys()
						if (t1 and t2):
							for species1 in taxalists[patient1][delimname][taxa]:
								for species2 in taxalists[patient2][delimname][taxa]:
									if species1 != species2:
										if i>0:
											#assess whether or not they are the same i-1 taxa 
											#meaning if family then make sure the genera are different
											lowertaxa1 = species1.split(delims[i-1])[1].split(';')[0]
											lowertaxa2 = species2.split(delims[i-1])[1].split(';')[0]
											if lowertaxa1 != lowertaxa2:
												connected = 0
												for arg in argdict.keys():
													if ((species1 in argdict[arg][patient1]) and (species2 in argdict[arg][patient2])):
														connected = 1
												if connected == 1:
													linkedcount = linkedcount + 1
												allcount = allcount + 1
										else: #this is when it is genera and you already know that everything should be compared
											connected = 0
											for arg in argdict.keys():
												if ((species1 in argdict[arg][patient1]) and (species2 in argdict[arg][patient2])):
													connected = 1
											if connected == 1:
												linkedcount = linkedcount + 1
											allcount = allcount + 1
					if patient1 == patient2:
						withinbetween = 'within'
					else:
						withinbetween = 'between'
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(patient1, patient2, linkedcount, allcount, withinbetween,delimname))


