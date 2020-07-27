
# WRK=/workdir/users/agk85/CDC2
# NAME=B314-2
# python contig_histogram_prep.py \
# -o $WRK/bins/${NAME}_mge_taxa_das_2_fgs \
# -i $WRK/bins/${NAME}_das_2_contigs_v_bins.txt \
# -c $WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
#files should come out as arg_orgcounts_${NAME}_species1.txt or species0.txt if not including hic

import argparse
from argparse import RawDescriptionHelpFormatter

#def getOptions():
description="""Programs reads through contig-bin files for a patient, pulls out the args and aggregates taxa at particular levels"""

parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-i','--infiles', dest="inhandles", action='store', required=True,  help='Infiles [REQUIRED]', metavar="INFILE")
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED]', metavar="OUTFILE")
parser.add_argument('-a','--mgeclust', dest="mgeclust", action='store', required=True, help='Arg file [REQUIRED]', metavar="ARGS")
parser.add_argument('-p','--patient', dest="patient", action='store', required=True, help='Patient [REQUIRED]', metavar="PATIENT")
parser.add_argument('-c','--contacts', dest="contacts", action='store', type=float, required=True, help='Min Contacts [REQUIRED]', metavar="CONTACTS")


args = parser.parse_args()


percent_thresh = 0.1
contact_thresh = args.contacts

levels = ['k__','p__','c__','o__','f__','g__', 's__']

mge_taxa_count = {}
#big dict
for level in levels:
	mge_taxa_count[level]={}


all_mgeids = []
with open(args.mgeclust) as mgefile:
	for line in mgefile:
		cluster = line.split('\t')[0]
		all_mgeids.append(cluster)



threshes = ['contacts','percentage']
bintypes = ['quality','anybin','bin+contig']
residencies = ['resident_hic','resident_nohic','hic']
#initiation
mge_taxa = {}
for thresh in threshes:
	mge_taxa[thresh] = {}
	for bintype in bintypes:
		mge_taxa[thresh][bintype]={}
		for residency in residencies:
			mge_taxa[thresh][bintype][residency] = {}
			for mgeid in all_mgeids:
				mge_taxa[thresh][bintype][residency][mgeid] = []
			

#iterate through file gathering bin taxa
#count   sample  contig  bin     bin_length      quality association_type        total_count     trans_count     norm_l  norm_rf contig_taxonomy kraken_taxonomy gtdb_taxonomy   arg_prsence     arg_clusters    arg_genes       mge_presence    mge_clusters    mge_genes
#cluster_resident, cluster_hic, contig_resident, contig_hic
#NA, BAD, LQ, MQ, HQ
mgelist = []
contact_count = 0
infiles = args.inhandles.split(',')
for inhandle in infiles:
	with open(inhandle) as infile:
		header = infile.readline()
		wonkyline = infile.readline()  ###############is this still a thing????
		for line in infile:
			mgepresence = line.split('\t')[17]
			if mgepresence == '1':
				mgeids = line.split('\t')[18]
				#temporary fix:
				if type(mgeids) is not list:
					mgeids = mgeids.split(',')
				for mgeid in mgeids:
					mgelist.append(mgeid)
					qualitybin = 0
					anybin = 1
					bin_or_contig = 0
					quality = line.split('\t')[5]
					try:
						prop = float(line.split('\t')[8])/float(line.split('\t')[7])
					except ZeroDivisionError:
						prop = 0
					taxonomy = line.split('\t')[12]
					contig_taxonomy = line.split('\t')[11]
					contig = line.split('\t')[2]
					contacts = float(line.split('\t')[8])
					association_type=line.split('\t')[6]
					#resident must be in the cluster or the contig of interest
					if association_type=='cluster_resident_withhicsupport':
						residency='resident_hic'
						if (quality == 'HQ' or quality == 'MQ'):
							mge_taxa['percentage']['quality'][residency][mgeid].append(taxonomy)
						if (quality == 'HQ' or quality == 'MQ' or quality == 'LQ'):
							mge_taxa['percentage']['anybin'][residency][mgeid].append(taxonomy)
						mge_taxa['percentage']['bin+contig'][residency][mgeid].append(taxonomy)				
						if (quality == 'HQ' or quality == 'MQ'):
							contact_count +=1
							mge_taxa['contacts']['quality'][residency][mgeid].append(taxonomy)
						if (quality == 'HQ' or quality == 'MQ' or quality == 'LQ'):
							mge_taxa['contacts']['anybin'][residency][mgeid].append(taxonomy)
						mge_taxa['contacts']['bin+contig'][residency][mgeid].append(taxonomy)		
					#if you are a resident, you get in for free	
					if association_type=='cluster_resident_withouthicsupport':
						residency='resident_nohic'
						if (quality == 'HQ' or quality == 'MQ'):
							mge_taxa['percentage']['quality'][residency][mgeid].append(taxonomy)
						if (quality == 'HQ' or quality == 'MQ'or quality == 'LQ'):
							mge_taxa['percentage']['anybin'][residency][mgeid].append(taxonomy)
						mge_taxa['percentage']['bin+contig'][residency][mgeid].append(taxonomy)
						if (quality == 'HQ' or quality == 'MQ'):
							mge_taxa['contacts']['quality'][residency][mgeid].append(taxonomy)
						if (quality == 'HQ' or quality == 'MQ' or quality == 'LQ'):
							mge_taxa['contacts']['anybin'][residency][mgeid].append(taxonomy)
						mge_taxa['contacts']['bin+contig'][residency][mgeid].append(taxonomy)
					#but if you care coming in through hic- you gotta have some more contact
					if association_type=='cluster_hic' :
						residency='hic'		 
						if prop >= percent_thresh:
							if (quality == 'HQ' or quality == 'MQ'):
								mge_taxa['percentage']['quality'][residency][mgeid].append(taxonomy)
							if (quality == 'HQ' or quality == 'MQ' or quality == 'LQ'):
								mge_taxa['percentage']['anybin'][residency][mgeid].append(taxonomy)
							mge_taxa['percentage']['bin+contig'][residency][mgeid].append(taxonomy)		
						if contacts >= contact_thresh:
							if (quality == 'HQ' or quality == 'MQ'):
								contact_count +=1
								mge_taxa['contacts']['quality'][residency][mgeid].append(taxonomy)
							if (quality == 'HQ' or quality == 'MQ' or quality == 'LQ'):
								mge_taxa['contacts']['anybin'][residency][mgeid].append(taxonomy)
							mge_taxa['contacts']['bin+contig'][residency][mgeid].append(taxonomy)
						

total_mges = len(set(mgelist))
mgeset = list(set(mgelist))

print(contact_count)
header = 'Count\tPatient\tThresh\tBintype\tResidency\tLevel\tARG\tTaxaCount\tTotal_mges\tTaxaset\tMin_contacts\n'
count = 0
with open(args.outhandle,'w') as outfile:
	outfile.write(header)
	for thresh in threshes:
		for bintype in bintypes:
			for level in levels:
				for residency in residencies:
					for mgeid in mgeset:
							taxalist = []
							for taxonomy in mge_taxa[thresh][bintype][residency][mgeid]:
								if taxonomy != 'NA' and taxonomy != '.':
									taxa = taxonomy.split(level)[1].split(';')[0]
									uptotaxa = taxonomy.split(level)[0] + level + taxonomy.split(level)[1].split(';')[0]
									if taxa != '':
										taxalist.append(uptotaxa)
							taxaset = set(taxalist)
							taxasetlen = len(taxaset)
							if taxasetlen > 0:
								outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(str(count), args.patient, thresh, bintype, residency, level, mgeid, str(taxasetlen), str(total_mges), ','.join(taxalist), str(args.contacts)))
								count += 1
