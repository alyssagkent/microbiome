#This file loops over files in bwa_alignments and grabs the rpkm if it is above 80% coverage spitting out a samp vs. arg table.
#python arg_sample_rpkm.py /workdir/users/agk85/CDC/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95 /workdir/users/agk85/CDC/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95.txt /workdir/users/agk85/CDC/arg_v_org/metagenomes3/mapping/bwa_alignments_95_95/arg_v_samp_95_95_names.txt

import glob
import collections
import sys

folder = sys.argv[1] #use the full path
outhandle = sys.argv[2]
outnamehandle = sys.argv[3]
filepaths = glob.glob(folder + '/*rpkm.txt')

#make a dictionary relating the protein with its antibiotic resistance gene name
arg_name_dict = collections.defaultdict(dict)
card_arg_name = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/card_prot_id_type.txt'
resfams_arg_name='/workdir/users/agk85/CDC/arg_v_org/metagenomes3/resfams_prot_id_type.txt'
with open(card_arg_name) as f:
	for line in f:
		r = line.split('\t')[0].strip()
		if r != 'ORF_ID':
			arg_name_dict[r]['name'] = line.strip().split('\t')[1]
			arg_name_dict[r]['card'] = line.strip().split('\t')[1]

with open(resfams_arg_name) as f:
	for line in f:
		r = line.split('\t')[0].strip()
		if r != 'ORF_ID':
			arg_name_dict[r]['name'] = line.strip().split('\t')[1]
			arg_name_dict[r]['resfams'] = line.strip().split('\t')[1]


#this assumes your rpkm files are all in the same order
#so if something funky happens there you will need to reassess

with open(outhandle,'w') as outfile:
	filehandle = filepaths[0]
	with open(outnamehandle,'w') as outnames:
		outnames.write('Protein\tName\tCARD\tResfams\n')
		with open(filehandle) as infile:
			header = infile.readline()
			for line in infile:
				arg = line.split(',')[0]
				if arg != '*':
					outfile.write(',' + arg)
					try:
						card = arg_name_dict[arg]['card']
					except:
						card = 'NA'
					try:
						resfams = arg_name_dict[arg]['resfams']
					except:
						resfams = 'NA'
					outnames.write('{0}\t{1}\t{2}\t{3}\n'.format(arg, arg_name_dict[arg]['name'] , card, resfams))
		outfile.write('\n')
	for filehandle in filepaths:
		with open(filehandle) as infile:
			header = infile.readline()
			sample = filehandle.split('/')[9].split('.')[0]
			outfile.write(sample)
			for line in infile:
				arg = line.split(',')[0]
				if arg != '*':
					rpkm = float(line.split(',')[1])
					coverage= float(line.strip().split(',')[2])
					if coverage >= 80:
						outfile.write(',' + str(rpkm))
					else:
						outfile.write(',0.0')
			outfile.write('\n')

	




