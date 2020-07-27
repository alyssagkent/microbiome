#this script figures how how many scaffolds have trans hic reads 
#either any hic-trans reads or only looking at pairs of scaffolds that have at least 2 readpairs
import glob

print('initializaiton')





samples = []
transcontigs2 = {}
transcontigs5 = {}
transcontigs = {}
validcontigs = {}
tbldict = {}
samplist = glob.glob('/workdir/users/agk85/CDC2/todo/*')
for samplong in samplist:
	sample = samplong.split('/')[-1]
	samples.append(sample)
	transcontigs2[sample] = []
	transcontigs5[sample] = []
	transcontigs[sample] = []
	validcontigs[sample] = []
	tbldict[sample] = []

samples.sort()
print('allvalidpairs')

#tbldict
tblpaths = glob.glob('/workdir/users/agk85/CDC2/combo_tables/metagenomes/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		sample = tblfile.split('/')[-1].split('_')[0]
		if sample in samples:
			header = t.readline()
			for line in t:
				tbldict[sample].append(line.split('\t')[0])


print('valid pairs')
infiles = glob.glob('/workdir/users/agk85/CDC2/hicpro/output/*_output/hic_results/data/*/*_allValidPairs')
for infile in infiles:
	print(infile)
        with open(infile) as f:
		for line in f:
			sample = infile.split('/')[-2]
			contig1 = line.split('\t')[1]
			contig2 = line.split('\t')[4]
			validcontigs[sample].append(contig1)
			validcontigs[sample].append(contig2)

print('trans pairs')
infiles = glob.glob('/workdir/users/agk85/CDC2/hicpro/output/*_output/hic_results/data/*/*_trans_primary_0_ncol_noeuks.txt')
for infile in infiles:
	print(infile)
	with open(infile) as f:
		for line in f:
			sample = infile.split('/')[-2]
			contig1 = line.split('\t')[0]
			contig2 = line.split('\t')[1]
			count = float(line.strip().split('\t')[2])
			transcontigs[sample].append(contig1)
			transcontigs[sample].append(contig2)
			if count>=2:
				transcontigs2[sample].append(contig1)
				transcontigs2[sample].append(contig2)
			if count>=5:
				transcontigs5[sample].append(contig1)
				transcontigs5[sample].append(contig2)

print('outputting')
with open("/workdir/users/agk85/CDC2/hicpro/figures/unique_scfs_hicmapped.txt",'w') as outfile:
	outfile.write("Sample\tTotalContigs\tValid1+\tTrans1+\tTrans2+\tTrans5+\n")
	for sample in samples:
		total = len(set(tbldict[sample]))
		validunique = len(set(validcontigs[sample]))
		tcunique = len(set(transcontigs[sample]))
		tc2unique =  len(set(transcontigs2[sample]))
		tc5unique = len(set(transcontigs5[sample]))
		outfile.write(sample + '\t' +str(total) + '\t' + str(validunique) + '\t' + str(tcunique) + '\t' + str(tc2unique) + '\t' + str(tc5unique)+ '\n')
