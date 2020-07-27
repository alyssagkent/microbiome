#remove_eukaryote_trans.py
import glob
import sys
euklist = []
folder = sys.argv[1]
pid = sys.argv[2]
name=sys.argv[3]
WRK='/workdir/users/agk85/' + folder
eukpaths = glob.glob(WRK + '/combo_tables/metagenomes/*_euk.txt')
for eukhandle in eukpaths:
	with open(eukhandle) as eukfile:
		for line in eukfile:
			euklist.append(line.split('\t')[0])

#changed from withexcise
transpaths = glob.glob(WRK + '/hicpro/output/*/hic_results/data/' + name + '/*_trans_primary_' + pid + '_ncol.txt')
for transpath in transpaths:
	transout  = transpath.split('.txt')[0] + '_noeuks' + '.txt'
	with open(transpath) as transfile:
		with open(transout,'w') as outfile:
			for line in transfile:
				r1 = line.split('\t')[0]
				r2 = line.split('\t')[1]
				if (r1 in euklist or r2 in euklist):
					continue
				else:
					outfile.write(line)


transpaths = glob.glob(WRK + '/hicpro/output/*/hic_results/data/' + name + '/*_trans_primary_' + pid + '_ncol_withexcise.txt')
for transpath in transpaths:
	transout  = transpath.split('.txt')[0] + '_noeuks' + '.txt'
	with open(transpath) as transfile:
		with open(transout,'w') as outfile:
			for line in transfile:
				r1 = line.split('\t')[0]
				r2 = line.split('\t')[1]
				if (r1 in euklist or r2 in euklist):
					continue
				else:
					outfile.write(line)


transpaths = glob.glob(WRK + '/hicpro/output/*/hic_results/data/' + name + '/*_trans_primary_' + pid + '.txt')
for transpath in transpaths:
	transout  = transpath.split('.txt')[0] + '_noeuks' + '.txt'
	with open(transpath) as transfile:
		with open(transout,'w') as outfile:
			for line in transfile:
				r1 = line.split('\t')[2]
				r2 = line.split('\t')[8]
				if (r1 in euklist or r2 in euklist):
					continue
				else:
					outfile.write(line)

#inhandle = WRK + '/hic/mapping/trans_primary_ncol_' + pid + '_withexcise.txt'
#outhandle = WRK+'/hic/mapping/trans_primary_ncol_' + pid + '_withexcise_noeuks.txt'
#deleted_trans = WRK+'/hic/mapping/trans_primary_ncol_' + pid + '_withexcise_euks.txt'



#inhandle = WRK + '/newhic/mapping/trans_primary_ncol_98.txt'
#outhandle = WRK+'/newhic/mapping/trans_primary_ncol_98_noeuks.txt'
#deleted_trans = WRK+'/newhic/mapping/trans_primary_ncol_98_euks.txt'
#with open(inhandle) as infile:
#	with open(outhandle,'w') as outfile:
#		with open(deleted_trans,'w') as deleted:
#			for line in infile:
#				r1 = line.split('\t')[0]
#				r2 = line.split('\t')[1]
#				if (r1 in euklist or r2 in euklist):
#					deleted.write(line)
#				else:
#					outfile.write(line)
				
