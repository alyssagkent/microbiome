import glob
from Bio import SeqIO
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
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
#good_patients = ['B370']

#tbldict
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC2/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()



bad_keys = []
#NB500947:511:HTL7JBGX3:1:11101:1732:1121        16      B314-1_scaffold_5744    18      0       81M     NB500947:511:HTL7JBGX3:1:11101:1732:1121        0       B314-1_scaffold_23951   1049    0       81M
infiles = glob.glob('/workdir/users/agk85/CDC2/todo/*')
#name = 'B335-1'
#label = 'B335-1'
for infile in infiles:
	mobile = {}
	linkcounter = {}
	cislinkcounter = {}
	name = infile.split('/')[-1]
	patient = name.split('-')[0]
	label = name
	if patient in good_patients:
		rpkm_dict = {}
		#normally I would do this but mapping was weird with the snakemake so I made a folder for them
		inrpkm = '/workdir/users/agk85/CDC2/mapping/' + name + '/' + name + '.rpkm.txt'
		inrpkm = '/workdir/users/agk85/CDC2/mapping/rpkms/' + name + '.rpkm.txt'
		with open(inrpkm) as rpkmfile:
			header = rpkmfile.readline()
			for line in rpkmfile:
				scfid = line.split(',')[0]
				rpkm_dict[scfid] = (line.split(',')[1], line.split(',')[2].strip())
		inscaffold= '/workdir/users/agk85/CDC2/prodigal_excise/metagenomes/' + name + '/' + name + '_scaffold.fasta'
		ref_dict = SeqIO.to_dict(SeqIO.parse(inscaffold, "fasta"))
		for scaffold in ref_dict.keys():
			linkcounter[scaffold] = 0
			cislinkcounter[scaffold]=0
		combotable = '/workdir/users/agk85/CDC2/combo_tables/metagenomes/'+name + '_master_scf_table.txt'

		with open(combotable) as ct:
			header = ct.readline()
			for line in ct:
				mobileness = line.split('\t')[-2]
				print(mobileness)
				scf = line.split('\t')[0]
				if mobileness == 'mge':
					mobile[scf] = 1
				else:
					mobile[scf] = 0
		#counts per link attached to each contig #change this to B314-1.nss.filter_98.combined.txt if you want all hits
		inhandle = '/workdir/users/agk85/CDC2/hicpro/output/'+ name +'_output/hic_results/data/'+name+'/'+ name +'_trans_primary_0_noeuks.txt'
		allhandle = '/workdir/users/agk85/CDC2/hicpro/output/'+ name +'_output/hic_results/data/'+name+'/'+ name +'_allValidPairs'
		with open(inhandle) as infile:
			for line in infile:
				contig1 = line.split('\t')[2]
				contig2 = line.split('\t')[8]
				linkcounter[contig1] += 1
				linkcounter[contig2] += 1
		with open(allhandle) as infile:
			for line in infile:
				contig1 = line.split('\t')[1]
				contig2 = line.split('\t')[4]
				if contig1 == contig2:
					cislinkcounter[contig1] +=1
	for key, h in linkcounter.iteritems():
		l = len(ref_dict[key])
		lengths.append(l)
		hic_hits.append(h)
		hic_cishits.append(cislinkcounter[key])
		try:
			mobility.append(mobile[key])
		except:
			mobility.append('0')
			bad_keys.append(key)
		names.append(name)
		scfids.append(key)
		rpkms.append(rpkm_dict[key][0])
		coverages.append(rpkm_dict[key][1])

outhandle = '/workdir/users/agk85/CDC2/read_distributions/All_length_hichits_mobility_abundance.txt'
header = 'Scfid\tName\tLength\tHic_hits\tMobility\tRPKM\tCoverage\tHic_cishits\n'
with open(outhandle, 'w') as outfile:
	outfile.write(header)
	for i in range(len(lengths)):
		n = names[i]
		l = lengths[i]
		h = hic_hits[i]
		a = hic_cishits[i]
		m = mobility[i]
		s = scfids[i]
		r = rpkms[i]
		c = coverages[i]
		#t = taxonomies[i]
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(s,n,l,h,m, r, c,a))

