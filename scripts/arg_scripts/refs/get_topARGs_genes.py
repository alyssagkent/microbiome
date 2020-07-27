#this script grabs top arg genes so the blast is way shorter than blasting EVERYTHING
from Bio import SeqIO
import glob
topargs = []
refhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/topargs/all_clusters_topargs_mergednames.txt'
with open(refhandle) as reffile:
	header = reffile.readline()
	for line in reffile:
		cluster = line.split('\t')[0]
		topargs.append(cluster)


samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

inhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/args_95_nr.fna.clstr'
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

recs =[]
seqhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/arg_prot.fasta'
for rec in SeqIO.parse(seqhandle,'fasta'):
	gene = rec.id
	cluster = protclusterdict[gene]
	if cluster in topargs:
		recs.append(rec)


SeqIO.write(recs, '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/patric_comparisons/toparg_prots.fasta','fasta')


	
