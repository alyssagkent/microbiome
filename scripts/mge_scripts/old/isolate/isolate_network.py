import igraph
from igraph import *
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt


isolate_associations = {}
clusters = {}
all_isolates = []
with open("mge_isolate_99_nr.fna.clstr.tbl") as infile:
	for line in infile:
		if 'NODE' in line.split('\t')[3]:
			prots = line.split('\t')[3].split(',')
			isolate_prots = []
			cluster = line.split('\t')[0]
			clusters[cluster] = []
			for prot in prots:
				if 'NODE' in prot:
					isolate_prots.append(prot)
			for isolateprot in isolate_prots:
				isolate = isolateprot.split('_')[0]
				all_isolates.append(isolate)
				clusters[cluster].append(isolate)
				for otherisolateprot in isolate_prots:
					otherisolate = otherisolateprot.split('_')[0]
					if isolate != otherisolate:
						try:
							isolate_associations[isolate +':' + otherisolate].append(cluster)
						except:
							isolate_associations[isolate +':' + otherisolate] = [cluster]


isolates = list(set(all_isolates))
with open("genes_vs_isolates.txt",'w') as outfile:
	for key in clusters:
		if len(clusters[key])>1:
			outfile.write('\t' + key)
	outfile.write('\n')
	for isolate in isolates:
		outfile.write(isolate)
		for cluster in clusters:
			if len(clusters[cluster])>1:
				if isolate in clusters[cluster]:
					outfile.write('\t1')
				else:
					outfile.write('\t0')
	 	outfile.write('\n')


outfilehandle = 'mge_isolate_99_ncol.txt'
with open(outfilehandle,'w') as outfile:
	for isolate1isolate2 in isolate_associations:
		clusters = isolate_associations[isolate1isolate2]
		isolate1 = isolate1isolate2.split(':')[0]
		isolate2 = isolate1isolate2.split(':')[1]
		count = str(len(clusters))
		outfile.write('{0}\t{1}\t{2}\n'.format(isolate1, isolate2, count))


