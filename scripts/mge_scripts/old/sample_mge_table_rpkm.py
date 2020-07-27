#This file loops over files in bwa_alignments and grabs the rpkm if it is above 80% coverage spitting out a samp vs. mge table.
#python mge_sample_rpkm.py /workdir/users/agk85/CDC/tables/metagenomes/mapping/bwa_alignments_mge_95_95 /workdir/users/agk85/CDC/tables/metagenomes/mapping/bwa_alignments_mge_95_95/mge_v_samp_95_95.txt 
import glob
import collections
import sys

folder = sys.argv[1] #use the full path
outhandle = sys.argv[2]
outnamehandle = sys.argv[3]
filepaths = glob.glob(folder + '/*rpkm.txt')

#this assumes your rpkm files are all in the same order
#so if something funky happens there you will need to reassess

with open(outhandle,'w') as outfile:
	filehandle = filepaths[0]
	with open(outnamehandle,'w') as outnames:
		outnames.write('Protein\n')
		with open(filehandle) as infile:
			header = infile.readline()
			for line in infile:
				clusterrep = line.split(',')[0]
				if clusterrep != '*':
					outfile.write(',' + clusterrep)
		outfile.write('\n')
	for filehandle in filepaths:
		with open(filehandle) as infile:
			header = infile.readline()
			sample = filehandle.split('/')[9].split('.')[0]
			outfile.write(sample)
			for line in infile:
				clusterrep = line.split(',')[0]
				if clusterrep != '*':
					rpkm = float(line.split(',')[1])
					coverage= float(line.strip().split(',')[2])
					if coverage >= 80:
						outfile.write(',' + str(rpkm))
					else:
						outfile.write(',0.0')
			outfile.write('\n')

	




