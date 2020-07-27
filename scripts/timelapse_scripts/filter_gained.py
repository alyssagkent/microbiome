import sys

#filter steps
#USAGE: python filter_gained.py timelapse_arg_org_2_alltaxalevels_gained.txt timelapse_arg_org_2_alltaxalevels_gained_filtered_bestconnections.txt timelapse_arg_org_2_alltaxalevels_gained_filteredindex.txt

goodlist ={}
with open(sys.argv[1]) as infile,open(sys.argv[2],'w') as outfile,open(sys.argv[4],'w') as outfile3:
	header = infile.readline()
	outfile.write(header)
	outfile3.write(header.strip() + '\tfilter\tfilter_reason\thigher_order\thigher_level_t1\thigher_level_t2\ttrans_taxa2_t1_morethan_0\ttrans_taxa2_t2_lessthan_2\tt1contigpres\ttrans__taxa1_t1_equals_0\tgenetype\n')
	for line in infile:
		count,connection_count,patient1,t1,t2,level,taxa1,taxa1_fulltaxa,ho,taxa1_t1_abund,taxa1_t2_abund,taxa2,fulltaxa2,taxa2_t1_abund,taxa2_t2_abund, gene,genename,t1rpkm,genecontigs_t2, t2rpkm,trans_taxa2_t1,trans_taxa2_t2,num_nont1,other_tps,contig_taxa_connection,contigfulltaxa1_t1,t1ll,trans_taxa1_t1,trans_taxa1_t2,num_bincontigs_t1, num_bincontigs_t2,binids_t1,binids_t2= line.strip().split('\t')
		if 'arg' in sys.argv[1]:
			genetype = 'ARG'
		else:
			genetype = 'MGE'	
		if (level == 's__' or level == 'g__') and ho=='0' and trans_taxa2_t1=='0' and (float(trans_taxa2_t2)>4) and (contig_taxa_connection == 't1contigabs') and (float(trans_taxa1_t1)>0) and (float(t1ll)>5):
			outfile3.write(line.strip()+'\tfiltered_good\tgood\t0\t0\t0\t0\t0\t0\t0\t' + genetype+'\n')
		else:
reason =[]
r level == 'g__') and ho=='0' and trans_taxa2_t1=='0' and (float(trans_taxa2_t2)>4) and (contig_taxa_connection == 't1contigabs') and (float(trans_taxa1_t1)>0) and (float(t1ll)>5):
			binary = []
			if ho != '0':
				reason.append('higher_order')
				binary.append('1')
			else:
				binary.append('0')
			if (float(t1ll)<5):
				reason.append('higher_level_t1')
				binary.append('1')
			else:
				binary.append('0')
			if (level!='s__' and level!='g__' and level!='f__'):
				reason.append('higher_level_t2')
				binary.append('1')
			else:
				binary.append('0')
			if (trans_taxa2_t1!='0'):
				reason.append('trans_taxa2_t1>0')
				binary.append('1')
			else:
				binary.append('0')
			if (float(trans_taxa2_t2)<2):
				reason.append('trans_taxa2_t2<2')
				binary.append('1')
			else:
				binary.append('0')
			if (contig_taxa_connection!= 't1contigabs'):
				reason.append('t1contigpres')
				binary.append('1')
			else:
				binary.append('0')
			if (float(trans_taxa1_t1)<1):
				reason.append('trans_taxa1_t1=0')
				binary.append('1')
			else:
				binary.append('0')
			outfile3.write(line.strip()+'\tfiltered_bad\t'+','.join(reason)+'\t' +'\t'.join(binary)+'\t' + genetype +'\n')
		if (level == 's__' or level == 'g__'):
			if ho == '0':
				if (trans_taxa2_t1 == '0' and float(trans_taxa2_t2)>4):
					if contig_taxa_connection == 't1contigabs':
						if(float(trans_taxa1_t1)>0):
							if(float(t1ll)>5):
								goodlist[connection_count] = line
								outfile.write(line)
		

with open(sys.argv[3],'w') as outfile2:
	outfile2.write(header)
	for connection_count in goodlist.keys():
		line = goodlist[connection_count]
		outfile2.write(line)
