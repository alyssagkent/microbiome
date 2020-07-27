#Goal per patient, per organism in that patient, assess how abundant they are from an average across all contigs with that annotation
#also across all time points

import glob
import collections
import sys
from best_org import best_org as bo

depth = '1'
#tbldict
tbldict = {}
tblpaths = glob.glob('/workdir/users/agk85/CDC/combo_tables/metagenomes4/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()


best_org_dict= {}
for node in tbldict:
	inline = tbldict[node]
	#get the last column which should always be the bestorg
	best_org = inline.split('\t')[-1]
	if 'k__Virus' in best_org:
		best_org = '.'
	best_org_dict[node] = best_org

#now replace things that have community memberships 
#change this if you move down to genus level!!!!!
if str(depth) != '0':
	networkpaths = glob.glob('/workdir/users/agk85/CDC' + '/networks/memberships/*_min' + str('2') + '_' + '98' + '_addon_species_membership.txt')
	for networkpath in networkpaths:
		print networkpath
		with open(networkpath) as network:
			for line in network:
				node = line.split('"')[1]
				species = line.split('"')[3]
				cluster = line.split('"')[4].split('\n')[0]
				#make sure you don't overwrite the names of the good things but this assumes the community names are going to be better, if you change how you create the communities, you might need to reassess this
				if species != '.' and species != 'NA':
					if ('k__Virus' not in species):
						best_org_dict[node] = species


#this is relinking phage
#best_org_dict = bo(tbldict, depth)

rpkmdict = {}
rpkmpaths = glob.glob('/workdir/users/agk85/CDC' + '/mapping/metagenomes3/bwa_alignments_scaffolds/*.rpkm.txt')
for rpkmhandle in rpkmpaths:
	with open(rpkmhandle) as rpkmfile:
		header = rpkmfile.readline()
		for line in rpkmfile:
			scfid = line.split(',')[0]
			rpkm = line.split(',')[1]
			covg = float(line.split(',')[2].strip())
			if covg >= 80:
				rpkmdict[scfid] = rpkm
			if covg <= 80:
				rpkmdict[scfid] = '0'


outhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/all_taxonomic_abundances_basehic_nophage.txt'
with open(outhandle,'w') as outfile:
	for node in best_org_dict.keys():
		if node != '':
			taxa= [best_org_dict[node]]
			mge = tbldict[node].split('\t')[-2]
			sample = node.split('_')[0]
			patient = sample.split('-')[0]
			rpkm = rpkmdict[node]
			for taxon in taxa:
				if taxon != '.':
					try:
						taxon = taxon.replace("'", ";")
					except AttributeError:
						print taxon
						print taxa 
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(node, sample, patient, rpkm, taxon))
