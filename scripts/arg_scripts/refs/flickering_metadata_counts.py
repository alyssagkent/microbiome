#bash scripts to get linecounts and trans reads
#wc -l */*_trans_primary_noeuks.txt > hic_trans_reads.txt
#for file in ~/wd/data/CDC/metagenomes/merged/unzip/*.1.fastq; do echo $file; total=$(wc -l $file); echo ${total} >> mgm_linecounts.txt; done
#for file in ~/wd/data/CDC/hic/merged/*.1.fastq; do echo $file; total=$(wc -l $file); echo ${total} >> hic_linecounts.txt; done


#then two scripts to get the assembly sizes before and after trimming---probably won't use the after
#to get the metagenome assembly sizes

# from Bio import SeqIO
# import glob
# paths = glob.glob('/workdir/users/agk85/CDC/idba_rerun/metagenomes/scaffolds/*.fasta')
# outhandle = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/mgm_assemblysize.txt"
# with open(outhandle, 'w') as outfile:
# 	for path in paths:
# 		sum = 0
# 		sample = path.split('/')[8].split('_scaffold')[0]
# 		print(sample)
# 		for rec in SeqIO.parse(path, "fasta"):
# 			sum = sum + len(rec)
# 		outfile.write(sample + '\t' + str(sum)+ '\n')			
# 
# 
# paths = glob.glob('/workdir/users/agk85/CDC/idba_rerun/metagenomes/*/scaffold.fa')
# outhandle = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/mgm_untrimmed_assemblysize.txt"
# with open(outhandle, 'w') as outfile:
# 	for path in paths:
# 		sum = 0
# 		sample = path.split('/')[7]
# 		print(sample)
# 		for rec in SeqIO.parse(path, "fasta"):
# 			sum = sum + len(rec)
# 		outfile.write(sample + '\t' + str(sum) + '\n')	
########################################################################################################

#merge 5 files re flickering
import glob
paths = glob.glob("/workdir/users/agk85/CDC/todo/*")
sampledict = {}
samples = []
for path in paths:
	sample = path.split('/')[6]
	sampledict[sample] = {}
	samples.append(sample)

f1 = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/mgm_linecounts.txt"
f2 = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/hic_linecounts.txt"
f3 = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/mgm_assemblysizes.txt"
f4 = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/mgm_untrimmed_assemblysizes.txt"
f5 = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/hic_trans_reads.txt"
o1 = "/workdir/users/agk85/CDC/arg_v_org/metagenomes3/flickering/metacounts_combined.txt"
header = "Sample\tPatient\tMgm_reads\tHic_reads\tMgm_assembly\tMgm_untrimmed_assmbly\tHic_trans_reads\n"

with open(f1) as infile:
	for line in infile:
		sample = line.split('unzip/')[1].split('.1.fastq')[0]
		count = float(line.split(' ')[0])/4.0
		sampledict[sample]['mgm_reads'] = str(count)

with open(f2) as infile:
	for line in infile:
		sample = line.split('merged/')[1].split('hic.1.fastq')[0]
		count = float(line.split(' ')[0])/4.0
		sampledict[sample]['hic_reads'] = str(count)

with open(f3) as infile:
	for line in infile:
		sample = line.split('\t')[0]
		count = line.strip().split('\t')[1]
		sampledict[sample]['mgm_assemb'] = count

with open(f4) as infile:
	for line in infile:
		sample = line.split('\t')[0]
		count = line.strip().split('\t')[1]
		sampledict[sample]['mgm_untrim_assemb'] = count

with open(f5) as infile:
	for line in infile:
		sample = line.strip().split(' ')[1].split('/')[0]
		count = line.strip().split(' ')[0]
		sampledict[sample]['trans_reads'] = count		

samples.sort()
with open(o1,'w') as outfile:
	outfile.write(header)
	for sample in samples:
		patient = sample.split('-')[0]
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(sample, patient, \
		sampledict[sample]['mgm_reads'],sampledict[sample]['hic_reads'], sampledict[sample]['mgm_assemb'], \
		sampledict[sample]['mgm_untrim_assemb'], sampledict[sample]['trans_reads']))

