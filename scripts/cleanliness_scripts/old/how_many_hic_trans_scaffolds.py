#this script figures how how many scaffolds have trans hic reads 
#either any hic-trans reads or only looking at pairs of scaffolds that have at least 2 readpairs
import glob
from igraph import *

samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

samples.sort()

scfs = []
min2scfs = []
g = Graph.Read_Ncol('/workdir/users/agk85/CDC/newhic/mapping/trans_primary_ncol_98.txt',weights=True, directed=True)
g.to_undirected(mode="collapse",combine_edges="sum")
scfs = g.vs['name']
#only keep the links greater than 2
g.es.select(weight_lt=2).delete()

g.vs.select(_degree_lt=1).delete()
min2scfs = g.vs['name']


#B370-3_3552|phage|1334|7624_83130
basescfs = []
for scf in scfs:
	if 'phage' in scf:
		basescf = scf.split('_')[0] + '_scaffold_' + scf.split('_')[1].split('|')[0]
		basescfs.append(basescf)
	else:
		basescfs.append(scf)

uniqbasescfs = list(set(basescfs))

basemin2scfs = []
for scf in min2scfs:
	if 'phage' in scf:
		basescf = scf.split('_')[0] + '_scaffold_' + scf.split('_')[1].split('|')[0]
		basemin2scfs.append(basescf)
	else:
		basemin2scfs.append(scf)

uniqbasemin2scfs = list(set(basemin2scfs))


with open("/workdir/users/agk85/CDC/newhic/mapping/unique_scfs_hicmapped.txt",'w') as outfile:
	outfile.write("Sample\tAnyhic\tMin2\n")
	for sample in samples:
		basecount = 0
		basemin2count = 0
		for scf in uniqbasescfs:
			if sample in scf:
				basecount += 1
		for scf in uniqbasemin2scfs:
			if sample in scf:
				basemin2count += 1
		outfile.write(sample + '\t' + str(basecount) + '\t' + str(basemin2count) + '\n')
