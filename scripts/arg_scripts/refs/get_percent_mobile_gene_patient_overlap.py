#percent mobile gene patient overlap?
#goal run through mobile genes, and figure out what % have more than one patient, more than 2 more than 3?
import glob
samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

samples.sort()

with open ('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/args_95_nr.fna.clstr.tbl') as infile:
	zero = 0
	mt0 = 0
	mt1 = 0
	mt2 = 0
	mt3 = 0
	mt4 = 0
	mt5 = 0
	mt6 = 0
	for line in infile:
		cluster = line.split('\t')[0]
		genes = line.split('\t')[2].split(',')
		pats = []
		for gene in genes:
			sample = gene.split('_')[0]
			if sample in samples:
				pats.append(gene.split('-')[0])
		patset = len(list(set(pats)))
		if patset == 0:
			zero += 1
		if patset>0:
			mt0 += 1
		if patset >1:
			mt1 += 1
		if patset >2:
			mt2 += 1
		if patset>3:
			mt3 += 1
		if patset>4:
			mt4 += 1
		if patset>5:
			mt5 += 1
		if patset>6:
			mt6 += 1


with open('/workdir/users/agk85/CDC/arg_v_org/metagenomes3/ARG_patient_overlap.txt','w') as outfile:
	outfile.write('Zero gene clusters:{0}\t{1}\n'.format(str(zero), str(zero)))
	outfile.write('1+ gene clusters:{0}\t{1}\n'.format(str(mt0), str(round(100*float(mt0)/mt0,4))))
	outfile.write('2+ gene clusters:{0}\t{1}\n'.format(str(mt1), str(round(100*float(mt1)/mt0,4))))
	outfile.write('3+ gene clusters:{0}\t{1}\n'.format(str(mt2), str(round(100*float(mt2)/mt0,4))))
	outfile.write('4+ gene clusters:{0}\t{1}\n'.format(str(mt3), str(round(100*float(mt3)/mt0,4))))
	outfile.write('5+ gene clusters:{0}\t{1}\n'.format(str(mt4), str(round(100*float(mt4)/mt0,4))))
	outfile.write('6+ gene clusters:{0}\t{1}\n'.format(str(mt5), str(round(100*float(mt5)/mt0,4))))
	outfile.write('7+ gene clusters:{0}\t{1}\n'.format(str(mt6), str(round(100*float(mt6)/mt0,4))))
	
	
			
