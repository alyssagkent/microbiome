import glob
from Bio import SeqIO
import numpy as np
from best_org_nophage import *

lengths = []
hic_hits = []
hic_cishits = []
mobility = []
names = []
scfids = []
rpkms = []
coverages = []

good_patients = ['B314','B316','B320','B331','B335','B357','B370','US3','US8']

#tbldict
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC2/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()

bincontigdict = {}
contigbindict = {}
with open('/workdir/users/agk85/CDC2/das/all_DASTool_scaffolds2bin.txt') as dasfile:
	for line in dasfile:
		contig = line.split('\t')[0]
		binid = line.strip().split('\t')[1]
		try:
			bincontigdict[binid].append(contig)
		except KeyError:
			bincontigdict[binid] = [contig]
		contigbindict[contig] = binid


#checkm quality
qualitydict = {}
cluster_length = {}
with open('/workdir/users/agk85/CDC2/das/all.stats') as checkm:
	for line in checkm:
		binid = line.split('\t')[0].split('.contigs')[0]
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
		qualitydict[binid] = quality
		cluster_length[binid] = length


transhiccontig = {}
validhiccontig = {}
for name in samples:
	inhandle = '/workdir/users/agk85/CDC2/hicpro/output/'+ name +'_output/hic_results/data/'+name+'/'+ name +'_trans_primary_0_noeuks.txt'
	allhandle = '/workdir/users/agk85/CDC2/hicpro/output/'+ name +'_output/hic_results/data/'+name+'/'+ name +'_allValidPairs'
	with open(inhandle) as infile:
		for line in infile:
			contig1 = line.split('\t')[2]
			contig2 = line.split('\t')[8]
			try:
				transhiccontig[contig1] += 1
			except KeyError:
				transhiccontig[contig1] = 1
			try:
				transhiccontig[contig2] += 1
			except KeyError:
				transhiccontig[contig2] += 1
	with open(allhandle) as infile:
		for line in infile:
			contig1 = line.split('\t')[1]
			contig2 = line.split('\t')[4]
			try:
				validhiccontig[contig1] +=1 
			except KeyError:
				validhiccontig[contig1] = 1
			try:
				validhiccontig[contig2] += 1
			except KeyError:
				validhiccontig[contig2] = 1


bincontigdict = {}
contigbindict = {}
transhicbinid = {}
validhicbinid = {}
transhictaxa = {}
validhictaxa = {}

with open('/wtrkdir/users/agk85/CDC2/das/all_DASTool_scaffolds2bin.txt') as dasfile:
	for line in dasfile:
		contig = line.split('\t')[0]
		binid = line.strip().split('\t')[1]
		try:
			transhicbinid[binid] += transhiccontig[contig]
		except KeyError:
			try:
				transhicbinid[binid] = transhiccontig[contig]
			except KeyError:
				transhicbinid[binid] = 0
		try:
                        validhicbinid[binid] += validhiccontig[contig]
                except KeyError:
                	validhicbinid[binid] = 0
		try:
			transhictaxa[taxa] += transhiccontig[contig]
		except KeyError:
			transhictaxa[taxa] = 0
		
			contigbindict[contig] = binid



outhandle = '/workdir/users/agk85/CDC2/read_distributions/All_bin_hic_counts.txt'
header = 'Binid\tHic_hits\tHic_cishits\tTaxonomy\n'
with open(outhandle, 'w') as outfile:
	outfile.write(header)
	for binid in binids:
		if qualitydict[binid] != 'BAD':
			trans = transhicbin[binid]
			valid = validhicbin[binid]
			kraken_taxa = krakenbin[binid]
			outfile.write('{0}\t{1}\t{2}\n'.format(binid,trans,valid,kraken_taxa))

