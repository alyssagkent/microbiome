#make the dictionary with how many reads are in the trimmed 99 one
import sys
name = sys.argv[1]
pid = sys.argv[2]
inhandle = '/workdir/users/agk85/CDC2/hicpro/output/'+name+'_output/bowtie_results/bwt2/'+name+'/'+name+'_'+name+'_scaffold.fasta.bwt2pairs_' + pid + '.sam'
vphandle = '/workdir/users/agk85/CDC2/hicpro/output/'+name+'_output/hic_results/data/'+name+'/'+name+'_allValidPairs'
outhandle = '/workdir/users/agk85/CDC2/hicpro/pid_counts/' + name + '_' + pid +'_pid_counts_withtrans.txt'

readdict = {}
with open(inhandle) as infile:
	for line in infile:
		readid = line.split('\t')[0]
		try:
			readdict[readid]+=1
		except:
			readdict[readid]=1

#open the file with all of the valid pairs
total = 0
goodcount = 0
transtotal = 0
goodtranscount = 0
with open(vphandle) as infile:
	for line in infile:
		total +=1
		readid = line.split('\t')[0]
		try:
			c = readdict[readid]
			if c == 2:
				goodcount +=1
			elif c != 1:
				print(c)
		except KeyError:
			a = 1
		contig1 = line.split('\t')[1]
		contig2 = line.split('\t')[4]
		if contig1 != contig2:
			transtotal +=1
			try:
				c = readdict[readid]
				if c == 2:
					goodtranscount +=1
				elif c != 1:
					print(c)
			except KeyError:
				a = 1

#go through each id, check if there are two for the above dictionary
#keep a count of all of them
#keep a count of the ones that pass
#print those that pass, those that are total, and the percentage in a oneline file
with open(outhandle,'w') as outfile:
	outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(name,str(goodcount),str(total),str(float(goodcount)/total),str(goodtranscount),str(transtotal),str(float(goodtranscount)/transtotal)))


