#"1","B357-2_scaffold_26863",110
#
#USAGE python split_cluster_hic.py /workdir/users/agk85/CDC/networks/memberships/file /workdir/users/agk85/CDC/networks/clusters/name name

import sys
inhandle = sys.argv[1]
reference = sys.argv[2]
outdir = sys.argv[3]
name = sys.argv[4]


from Bio import SeqIO

ref_dict = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))

clusterdict = {}
with open(inhandle) as infile:
	header = infile.readline()
	for line in infile:
		fields = line.strip().split('"')
		scf = fields[1]
		group = fields[4].split(',')[1]
		try:
			clusterdict[group].append(scf)
		except KeyError:
			clusterdict[group]=[scf]


for cluster in clusterdict.keys():
	scfs = clusterdict[cluster]
	scfrecs = []
	for scf in scfs:
		scfrecs.append(ref_dict[scf])
	
	SeqIO.write(scfrecs, outdir + '/Cluster_' + cluster + '.fasta','fasta')
			



