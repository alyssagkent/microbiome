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

taxataxa_dict = {}
for patient in patients:
	taxataxa_dict[patient] = {}
	for level in levels:
		taxataxa_dict[patient][level]= {}


taxa_color = {}
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
				level = fields[5]
				taxa = fields[9].split(',')
				taxaset = set(taxa)
				for taxon1 in taxaset:
					try:
						phylum = taxon1.split('p__')[1].split(';')[0]
					except:
						phylum = ''
					try:
						t = taxa_color[phylum]
					except KeyError:
						taxa_color[phylum] = '#%02X%02X%02X' % (r(),r(),r())
					for taxon2 in taxaset:
						if taxon1!=taxon2:
							try:
								taxataxa_dict[patient][level][taxon1+'|'+taxon2] += 1
								taxataxa_dict[patient][level][taxon2+'|'+taxon1] += 1
							except KeyError:
								taxataxa_dict[patient][level][taxon1+'|'+taxon2] = 1
								taxataxa_dict[patient][level][taxon2+'|'+taxon1] = 1


#output: Count	Patient	Level	taxa1	taxa2	ARGs	MGEs
header = 'Count\tPatient\tLevel\tTaxa1\tTaxa2\tGenecount\tBase1\tBase2\tColor1\tColor2\n'
count = 0
with open(args.outhandle,'w') as outfile:
	outfile.write(header)
	for patient in patients:
		for i in range(len(levels)):
			level = levels[i]
			levelname = levels[i]
			for key in taxataxa_dict[patient][level]:
				taxa1 = key.split('|')[0]
				taxa2 = key.split('|')[1]
				genecount = str(taxataxa_dict[patient][level][key])
				t1 = taxa1.split('; ')[-1].split(';')[0]
				t2 = taxa2.split('; ')[-1].split(';')[0]
				try:
					phylum1 = taxa1.split('p__')[1].split(';')[0]
				except:
					phylum1 = ''
				try:
					phylum2 = taxa2.split('p__')[1].split(';')[0]
				except:
					phylum2 = ''				
				c1 = '"' + taxa_color[phylum1] + '"'
				c2 = '"' + taxa_color[phylum2] + '"'
				outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(str(count), patient, level, taxa1, taxa2, genecount,t1,t2,c1,c2))
				count +=1


		