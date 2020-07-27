#make a more parse-able table from cd-hit output
#real specific right now for the mge's
import sys
inhandle = sys.argv[1]
#inhandle = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/mge_' + '99' + '_nr.fna.clstr'
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
			try:
					clusterprotdict[cluster].append(prot)
			except:
					clusterprotdict[cluster] = [prot]
			if '*' in line:
				clusternum_map[cluster] = prot


#make a tab delimited file in this format
#cluster number\tstarredid\tid1,id2,id3,id4\tnumber_of_proteins\tnumber_of_unique_patients\n
outhandle = inhandle + '.pattbl'
with open(outhandle, 'w') as outfile:
	for key in clusterprotdict.keys():
		prots = clusterprotdict[key]
		starred = clusternum_map[key]
		patients = set([prot.split('-')[0] for prot in prots])
		line = '{0}\t{1}\t{2}\t{3}\t{4}\n'.format(key, starred, ','.join(prots), str(len(prots)),str(len(patients)))
		outfile.write(line)
