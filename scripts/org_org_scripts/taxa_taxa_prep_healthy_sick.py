import argparse
from argparse import RawDescriptionHelpFormatter
from copy import deepcopy
import random

#def getOptions():
description="""Program turns contig_v_bin into taxa-taxa ARG and MGE counts"""

parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i','--in', dest="inhandle", action='store', required=True, help='Inhandle [REQUIRED]', metavar="MIN_PROP")
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED]', metavar="OUTFILE")
#parser.add_argument('-g','--gtdb', dest="gtdb_clusters", action='store', required=True, help='GTDB cluster taxonomies [REQUIRED]', metavar="GTDB")

args = parser.parse_args()
levels = ['k__','p__','c__','o__','f__','g__', 's__']
levelnames = ['kingdom','phylum','class','order','family','genus','species']
patients = ['B314','B316','B320','B331','B335','B357','B370','US3','US8','all']
h = ['US3','US8']
s = ['B314','B316','B320','B331','B335','B357','B370']


taxataxa_dict = {}
for patient in patients:
	taxataxa_dict[patient] = {}
	for level in levels:
		taxataxa_dict[patient][level]= {}

#init healthysickdict
hs_dict = {}
statuses = ['healthy','sick']
for status in statuses:
	hs_dict[status] = {}
	for level in levels:
		hs_dict[status][level] = {}

#way to get random colors
r = lambda: random.randint(0,255)

#Together_das_5_argtaxa_together.txt
#Count   Patient Thresh  Bintype Residency       Level   ARG     TaxaCount       Total_args      Taxaset Min_contacts
with open(args.inhandle) as infile:
	header = infile.readline()
	for line in infile:
		fields = line.strip().split('\t')
		thresh = fields[2]
		bintype = fields[3]
		if (thresh == 'contacts' and bintype == 'anybin'): #assuming this is what you want
			taxacount = int(fields[7])
			if taxacount>1:
				patient = fields[1]
				if patient != 'all':
					if patient in h:
						status = 'healthy'
					if patient in s:
						status = 'sick'
					level = fields[5]
					gene = fields[6]
					taxa = fields[9].split(',')
					taxaset = set(taxa)
					for taxon1 in taxaset:
						for taxon2 in taxaset:
							if taxon1!=taxon2:
								try:
									hs_dict[status][level][taxon1+'|'+taxon2].append(gene)
									hs_dict[status][level][taxon2+'|'+taxon1].append(gene)
								except KeyError:
									hs_dict[status][level][taxon1+'|'+taxon2]=[gene]
									hs_dict[status][level][taxon2+'|'+taxon1]=[gene]
								phyla1 = taxon1.split('; c__')[0]
								phyla2 = taxon2.split('; c__')[0]
								if phyla1 == phyla2:
									try:
										hs_dict[status]['p__'][phyla1+'|'+phyla2].append(gene)
									except KeyError:
										hs_dict[status]['p__'][phyla1+'|'+phyla2] =[gene]
									
header = 'Count\tPatient\tLevel\tTaxa1\tTaxa2\tGenecount\tBase1\tBase2\tGenes\n'
count = 0
with open(args.outhandle,'w') as outfile:
	outfile.write(header)
	for status in statuses:
		for i in range(len(levels)):
			finished =[]
			level = levels[i]
			levelname = levels[i]
			for orgorg in hs_dict[status][level].keys():
				taxa1 = orgorg.split('|')[0]
				taxa2 = orgorg.split('|')[1]
				if taxa1+taxa2 not in finished:
					finished.append(taxa2+taxa1)
					genecount = len(set(hs_dict[status][level][orgorg]))
					genes = ','.join(set(hs_dict[status][level][orgorg]))
					t1 = taxa1.split('; ')[-1].split(';')[0]
					t2 = taxa2.split('; ')[-1].split(';')[0]
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(str(count), status, level, taxa1, taxa2, genecount,t1,t2,genes))
					count +=1
