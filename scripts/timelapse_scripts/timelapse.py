import numpy as np
import glob
import collections
from Bio import SeqIO
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
 	"""Get arguments"""
 	description="""This script can be used to get arg versus organism capturing orgs on the 
 	contig and contigs up to N links away via Hi-C reads"""
 	parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
	parser.add_argument('-gc','--genecluster',dest='genecluster',action='store',required=True,type=str, help='Genes clstr file [Required]', metavar='ARGPID')
	parser.add_argument('-g','--genetype',dest='genetype',action='store',required=True,type=str, help='Gene [Required]', metavar='GENE')
	parser.add_argument('-m', '--minreads', dest='minreads', action='store', required=True, type=int, help='Min reads [Required]')
	parser.add_argument('-c','--connections', dest="connections", action='store', required=True, help='Contig vs bin file [REQUIRED]', metavar="CONNECTIONS")
	parser.add_argument('-i','--in', dest="cvb", action='store', required=True, help='Contig vs bin file [REQUIRED]', metavar="OUTFILE")
	parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="INFILE")
	args = parser.parse_args()
 	return(args)

args = getOptions()
####################################
#ok simpler you just want to know for each patient, out of all times when there is a connection, how consistent is it across the timepoints

####BIN STUFF######

levels = ['k__','p__','c__','o__','f__','g__', 's__']
samplelist = []
orgdict = {}
taxons = []
with open(args.binhandle) as binfile:
	for line in binfile:
		binid,patient,sample,quality,taxonomy = line.strip().split('\t')
		samplelist.append(sample)
		try:
			a = orgdict[sample]
		except KeyError:
			orgdict[sample] = {}
		if taxonomy != '.':
			for level in levels:
				taxa = taxonomy.split(level)[1].split(';')[0]
			if taxonomy != '.':
				for level in levels:
					taxa = taxonomy.split(level)[1].split(';')[0]
					orgdict[sample][taxa] = 1
					taxons.append(taxa)


samples = list(set(samplelist))
samples.sort()


print('Initialize the genefulldict')
genedict = {}
for sample in samples:
 	genedict[sample] = {}

with open(args.genecluster) as genefile:
	for line in genefile:
		cluster = line.split('\t')[0]
		genelist = line.split('\t')[2].split(',')
		for gene in genelist:
			sample = gene.split('_')[0]
			try:
				genedict[sample][cluster] = 1
			except KeyError:
				genedict[sample] = {}
				genedict[sample][cluster] = 1	



print('Initialize connectiondict')
connection_dict = {}
for sample in samples:
	connection_dict[sample] = {}
	for level in levels:
		connection_dict[sample][level] = []

counter = 0
#B314    B314-1  arg     2       5329    k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia_coli;
with open(args.connections) as infile:
	header = infile.readline()
	for line in infile:
		counter += 1
		patient,sample,genetype,minreads,geneid,taxonomy = line.strip().split('\t')
		#only if it has a species
		for level in levels:
			taxa = taxonomy.split(level)[1].split(';')[0]
			if taxa != '':
				connection_dict[sample][level].append((geneid,taxa,taxonomy,counter))

#get the arg rpkm
rpkmdict = {}
generepdict = {}
genenamedict = {}
if genetype == 'arg':
	with open('/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech_resfinder.txt') as infile:
		for line in infile:
			cluster = line.split('\t')[0]
			repgene = line.split('\t')[1]
			name = line.split('\t')[7] #give it the resfinder, then card then resfams
			if name == 'NA':
				name = line.split('\t')[3] #card
			if name == 'NA':
				name = line.split('\t')[2] #should be resfams
			genenamedict[cluster] = name
			generepdict[repgene] = cluster
	
	with open('/workdir/users/agk85/CDC2/args/mapping/bwa_alignments_99_99/arg_v_samp_99_99.txt') as rpkmfile:
		header = rpkmfile.readline()
		geneids = header.strip().split(',')
		sid = geneids.pop(0)
		for line in rpkmfile:
			rpkms = line.strip().split(',')
			sample = rpkms.pop(0)
			rpkmdict[sample] = {}
			for repgene, rpkm in zip(geneids, rpkms):
				clusterid = generepdict[repgene]
				rpkmdict[sample][clusterid] = rpkm
if genetype == 'mge':
	with open('/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.fna.clstr.tbl.annot') as infile:
		for line in infile:
			cluster = line.split('\t')[0]
			repgene = line.split('\t')[12]
			pfam = line.split('\t')[1]
			name = line.split('\t')[3] + '_'+ pfam #this is the MGE_combo_type + '_' + pfam annotation
			genenamedict[cluster] = name
			generepdict[repgene] = cluster
	with open('/workdir/users/agk85/CDC2/mobile/metagenomes/mapping/bwa_alignments_mge_99_99/mge_v_samp_99_99.txt') as rpkmfile:
		header = rpkmfile.readline()
		geneids = header.strip().split(',')
		sid = geneids.pop(0)
		for line in rpkmfile:
			rpkms = line.strip().split(',')
			sample = rpkms.pop(0)
			rpkmdict[sample] = {}
			for repgene, rpkm in zip(geneids, rpkms):
				clusterid = generepdict[repgene]
				rpkmdict[sample][clusterid] = rpkm

#get the taxa abundance
abunddict = {}
for sample in samples:
	abunddict[sample] = {}
with open('/workdir/users/agk85/CDC2/das/all_bintables_metaphlan.txt') as infile:
	for line in infile:
		binid,patient,sample,quality,taxonomy,k,p,c,o,f,g,s,best_binid_level,best_abund_level,diff = line.strip().split('\t')
		taxabunds = [k,p,c,o,f,g,s]
		if taxonomy == '.':
			taxonomy = 'k__; p__; c__; o__; f__; g__; s__;'
		for i in range(len(levels)):
			abund = taxabunds[i]
			level = levels[i]
			taxon = taxonomy.split(level)[1].split(';')[0]
			if taxon != '':
				abunddict[sample][taxon] = abund


#define the singleton connections:
goodbin_anyreads ={}
goodbin_anyreads_fulltaxa = {}
gene_taxa_trans = {}
contig_anyreads = {}
contig_anyreads_fulltaxa={}
gene_contig = {}
gene_contig_notaxon = {}
num_bincontigs = {}
gene_binid = {}
for sample in samples:
	goodbin_anyreads[sample] = {}
	goodbin_anyreads_fulltaxa[sample] = {}
	contig_anyreads[sample] = {}
	contig_anyreads_fulltaxa[sample]={}
	gene_taxa_trans[sample] = {}
	gene_contig[sample] = {}
	gene_contig_notaxon[sample] = {}
	num_bincontigs[sample] = {}
	gene_binid[sample] = {}
	for level in levels:
		goodbin_anyreads[sample][level] = {}
		goodbin_anyreads_fulltaxa[sample][level] = {}
		contig_anyreads[sample][level] = {}
		contig_anyreads_fulltaxa[sample][level] = {}
		gene_taxa_trans[sample][level] = {} 
		gene_contig[sample][level] = {}
		gene_contig_notaxon[sample][level] = {}
		num_bincontigs[sample][level] = {}
		gene_binid[sample][level] = {}

with open('/workdir/users/agk85/CDC2/bins/all_das_1_contigs_v_bins_contignum.txt') as dasfile:
	header = dasfile.readline()
	for line in dasfile:
		count,sample,contig,binid,bin_length,quality,association_type,total_count,transcount,norm_l,norm_rf,contig_taxonomy,kraken_taxonomy,gtdb_taxonomy,arg_presence,arg_clusters,arg_genes,mge_presence,mge_clusters,mge_genes,contig_num=line.strip().split('\t')
		trans_count = float(transcount)
		for level in levels:
			if kraken_taxonomy != 'NA': #this means that it is not an unannotated contig
				taxon = kraken_taxonomy.split(level)[1].split(';')[0]
				if genetype == 'arg':
					geneclusters = arg_clusters.split(',')
				if genetype == 'mge':
					geneclusters = mge_clusters.split(',')
				for cluster in geneclusters:
					try:
						num_bincontigs[sample][level][cluster][taxon]+=float(contig_num)
					except KeyError:
						try:
							num_bincontigs[sample][level][cluster][taxon]=float(contig_num)
						except KeyError:
							num_bincontigs[sample][level][cluster]={taxon:float(contig_num)}
					if 'contig' in association_type:
						try:
							contig_anyreads[sample][level][cluster].append(taxon)
						except KeyError:
							contig_anyreads[sample][level][cluster] = [taxon]
						try:
							contig_anyreads_fulltaxa[sample][level][cluster].append(kraken_taxonomy)
						except KeyError:
							contig_anyreads_fulltaxa[sample][level][cluster]=[kraken_taxonomy]
					if quality != 'NA':
						try:
							goodbin_anyreads[sample][level][cluster].append(taxon)
						except KeyError:
							goodbin_anyreads[sample][level][cluster] = [taxon]
						try:
							goodbin_anyreads_fulltaxa[sample][level][cluster].append(kraken_taxonomy)
						except KeyError:
							goodbin_anyreads_fulltaxa[sample][level][cluster] = [kraken_taxonomy]
						#count the trans reads
						try:
							gene_taxa_trans[sample][level][cluster][taxon]+=trans_count
						except KeyError:
							try:
								gene_taxa_trans[sample][level][cluster][taxon]=trans_count
							except KeyError:
								gene_taxa_trans[sample][level][cluster] = {taxon:trans_count}
						#add the contig that carries the cluster
						try:
							gene_contig[sample][level][cluster][taxon].append(contig)
						except KeyError:
							try:
								gene_contig[sample][level][cluster][taxon]=[contig]
							except KeyError:
								gene_contig[sample][level][cluster] = {taxon:[contig]}
						try:
							gene_contig_notaxon[sample][level][cluster].append(contig)
						except KeyError:
							gene_contig_notaxon[sample][level][cluster]=[contig]
						try:
							gene_binid[sample][level][cluster][taxon].append(binid)
						except KeyError:
							try:
								gene_binid[sample][level][cluster][taxon]=[binid]
							except KeyError:
								gene_binid[sample][level][cluster] = {taxon:[binid]}


#print(num_bincontigs)
timepoint_dict = {"B314-1":"1","B314-2":"2","B314-3":"3","B314-4":"4","B316-1":"1","B316-2":"2","B316-3":"3","B316-4":"4","B316-5":"5","B316-6":"6","B316-7":"7","B320-1":"1","B320-2":"2","B320-3":"3","B320-5":"4","B331-1":"1","B331-2":"2","B331-3":"3","B331-4":"4","B335-1":"1","B335-2":"2","B335-3":"3","B357-1":"1","B357-2":"2","B357-3":"3","B357-4":"4","B357-5":"5","B357-6":"6","B370-1":"1","B370-2":"2","B370-3":"3","B370-4":"4","B370-5":"5","US3-8":"1","US3-10":"2","US3-12":"3","US3-14":"4","US3-16":"5","US8-1":"1","US8-2":"2","US8-3":"3","US8-4":"4","US8-5":"5"}

def share_gene(t1,t2,connection):
	geneid = connection[0]
	if (geneid in genedict[t1] and geneid in genedict[t2]):
		sharedgene=1
	else:
		sharedgene=0
	return sharedgene

def share_taxa(t1,t2,connection):
	taxa = connection[1]
	if (taxa in orgdict[t1] and taxa in orgdict[t2]):
		sharedtaxa=1
	else:
		sharedtaxa=0
	return sharedtaxa


def get_connection_stats(t1,t2,connect1, connect2):
	shared = len(connect1 & connect2)
	uniq1 = connect1 - connect2 #those unique to t1
	uniq2 = connect2 - connect1 #those unique to t2
	unconnected1 = 0
	unconnected2 = 0
	uniqelement1 = 0
	uniqelement2 = 0
	gainlist = []
	lostlist = []
	gainsources = []
	#go thru elements unique to t1 (why unique?)
	for connection in uniq1:
		if share_gene(t1,t2,connection) and share_taxa(t1,t2,connection):
			unconnected1 += 1
			lostlist.append(connection)
		else:
			uniqelement1 += 1
	#same for t2
	for connection in uniq2:
		if share_gene(t1,t2,connection) and share_taxa(t1,t2,connection):
			unconnected2 += 1
			gainlist.append(connection)
			#check to see if the gene had any connection in the orig timepoint
			for con in connection_dict[t1]:
				if con[0] == connection[0]:
					gainsources.append((con+(connection[1],'dummy')))
		else:
			uniqelement2 += 1		
	tup =(shared,uniqelement1,uniqelement2,unconnected1,unconnected2,gainlist,lostlist,gainsources)
	return tup

# goal you need to get fulltaxa2 if it is 7 then work your way up from there
# the problem is that you have to fill in the rest with empties
emptylist = ['; p__; c__; o__; f__; g__; s__;', '; c__; o__; f__; g__; s__;','; o__; f__; g__; s__;','; f__; g__; s__;','; g__; s__;','; s__;',';']
#new approach, compare each full taxa with fulltaxa2 level by level, if there is no conflict, then that's a higher_order raised
def has_higher_order(fulltaxa2,fulltaxat1list): #I"m assuming that taxa2 !=''
	higher_order = 0
	for taxonomy1 in fulltaxat1list:
		taxa_conflict = 0
		for level in levels:
			taxa1 = taxonomy1.split(level)[1].split(';')[0]
			taxa2 = fulltaxa2.split(level)[1].split(';')[0]
			if (taxa1!=taxa2 and taxa1!= '' and taxa2 != ''):
				taxa_conflict = 1
		if taxa_conflict == 0:
			higher_order =1
	return higher_order

def lowest_level(taxonomy_list):
	lowest_level = -1 #empty
	for taxonomy in taxonomy_list:
		for i in range(len(levels)):
			level=levels[i]
			taxon = taxonomy.split(level)[1].split(';')[0]
			if taxon != '' and lowest_level < i:
				lowest_level = i
	return lowest_level+1

level='species'
count = 0
print('Running through all the combinations')
outhandle = '/workdir/users/agk85/CDC2/bins/timelapse/timelapse_{0}_org_{1}_alltaxalevels.txt'.format(args.genetype, level, str(args.minreads))
outhandle2 = '/workdir/users/agk85/CDC2/bins/timelapse/timelapse_{0}_org_{1}_alltaxalevels_gained.txt'.format(args.genetype, str(args.minreads))
with open(outhandle, 'w') as outfile, open(outhandle2,'w') as outfile2:
	header = 'T1\tT2\tLevel\tTimepoint_of_comparison\tTotal_T1\tTotal_T2\tShared\tUniqelement1\tUniqelement2\tUnconnected1\tUnconnected2\n'
	headerfull = 'Count\tConnection_count\tPatient\tT1\tT2\tLevel\tTaxa1\tFulltaxa1\tTaxa2_in_higher_order_taxa1s\tTaxa1_t1_abund\tTaxa1_t2_abund\tTaxa2\tFulltaxa2\tTaxa2_t1_abund\tTaxa2_t2_abund\tGeneid\tGene_name\tGene_contig\tGene_t1_abund\tGene_t2_abund\tTrans_taxa2_t1\tTrans_taxa2_t2\tNumber_non_t1_tps\tOther_tps\tContig_taxa_connection\tContig_taxa\tLowest_level_t1\tTrans_taxa1_t1\tTrans_taxa1_t2\tNum_bincontigs_t1\tNum_bincontigs_t2\tBinid_t1\tBinid_t2\n'
	outfile.write(header)
	outfile2.write(headerfull)
	for levelindex in range(len(levels)):
		level = levels[levelindex]
		print(level)
		for t1 in samples:
			patient1 = t1.split('-')[0]
			for t2 in samples:
				patient2 = t2.split('-')[0]
				if patient1 == patient2:
					connect1 = set(connection_dict[t1][level])
					connect2 = set(connection_dict[t2][level])
					tott1 = len(connect1)
					tott2 = len(connect2)
					shared,uniqelement1, uniqelement2,unconnected1,unconnected2,gainlist,lostlist,gainsources=get_connection_stats(t1,t2,connect1, connect2)
					outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format(patient1,t2,level,timepoint_dict[t1],tott1,tott2,shared,uniqelement1, uniqelement2,unconnected1,unconnected2))
					if timepoint_dict[t1] == '1':
						for connection in gainlist:
							connection_count=str(connection[3])
							gene = connection[0]
							taxa = connection[1]
							t1clusterpres = 't1clusterabs'
							t1contigpres = 't1contigabs'
							t1clustertaxa = []
							t1clusterfulltaxa = []
							taxa1_t1_abunds = []
							taxa1_t2_abunds = []
							taxa2_t1_abunds = []
							taxa2_t2_abunds = []
							t1contigtaxa = []
							t1contigfulltaxa = []
							try:
								if taxa in goodbin_anyreads[t1][level][gene]:
									t1clusterpres = 't1clusterpres'
								t1clustertaxa = t1clustertaxa + goodbin_anyreads[t1][level][gene] # [t.split('s__')[1] for t in goodbin_anyreads[t1][gene]]
								t1clusterfulltaxa = t1clusterfulltaxa + goodbin_anyreads_fulltaxa[t1][level][gene]
							except KeyError:	
								a = 1
							if taxa not in t1clustertaxa:
								try:
									if taxa in contig_anyreads[t1][level][gene]:
										t1contigpres = 't1contigpres'
									t1contigtaxa = t1contigtaxa + contig_anyreads[t1][level][gene] # [t.split('s__')[1] for t in contig_anyreads[t1][gene] if t !='NA']
									t1contigfulltaxa = t1contigfulltaxa + contig_anyreads_fulltaxa[t1][level][gene]
								except KeyError:
									a = 1
								number_nont1 = 1
								nont1 = [t2]
								for t3 in samples:
									patient3 = t3.split('-')[0]
									if (t3 != t1) and (t3 not in nont1) and (patient3 != patient2):
										if connection in connection_dict[t3]:
											number_nont1+=1
											nont1.append(t3)
								#Get gene rpkms
								t1rpkm = str(rpkmdict[t1][gene])
								t2rpkm = str(rpkmdict[t2][gene])
								#get the abundances
								taxa1_t1_abunds = []
								taxa1_t2_abunds = []
								for t1taxon in t1clustertaxa:
									try:
										taxa1_t1_abunds.append(abunddict[t1][t1taxon])
									except KeyError:
										taxa1_t1_abunds.append('NA')
									try:
										taxa1_t2_abunds.append(abunddict[t2][t1taxon])
									except:
										taxa1_t2_abunds.append('NA')
								taxa2_t1_abund = abunddict[t1][taxa]
								taxa2_t2_abund = abunddict[t2][taxa]
								taxa1 = ','.join(t1clustertaxa) #get the taxa in a csv list
								taxa1_fulltaxa = ','.join(t1clusterfulltaxa) #get the full taxa
								taxa1_t1_abund = ','.join(taxa1_t1_abunds)
								taxa1_t2_abund = ','.join(taxa1_t2_abunds)
								if taxa1 == '':
									taxa1 = 'NA_or_unannotated'
									taxa1_t1_abund = 'NA'
									taxa1_t2_abund = 'NA'
								if taxa1_fulltaxa=='':
									taxa1_fulltaxa='No_bin?'
								taxa2 = taxa
								fulltaxa2 = connection[2]
								genename = genenamedict[gene]
								try:
									trans_taxa2_t1 = gene_taxa_trans[t1][level][gene][taxa2]
								except KeyError:
									trans_taxa2_t1 = 0
								trans_taxa2_t2 = gene_taxa_trans[t2][level][gene][taxa2]
								try:
									trans_taxa1_t1 = gene_taxa_trans[t1][level][gene][taxa1]
								except KeyError:
									trans_taxa1_t1 = 0
								try:
									trans_taxa1_t2 = gene_taxa_trans[t2][level][gene][taxa1]
								except KeyError:
									trans_taxa1_t2 = 0
								num_nont1 = number_nont1
								other_tps = ','.join(nont1)
								#allgenecontigs_t1 = set(gene_contig_notaxon[t1][level][gene])
								#allgenecontigs_t2=set(gene_contig_notaxon[t2][level][gene])
								#allgenecontigs_taxonomies_t1 = ','.join([kraken_contig_taxonomy[scf] for scf in allgenecontigs_t1])
								#allgenecontigs_taxonomies_t2 = ','.join([kraken_contig_taxonomy[scf] for scf in allgenecontigs_t2])
								genecontigs_t2 = ','.join(set(gene_contig[t2][level][gene][taxa2]))
								genebinids_t2 = ','.join(set(gene_binid[t2][level][gene][taxa2]))
								temp = []
								for t1_taxon in t1clustertaxa:
									temp.append(','.join(gene_binid[t1][level][gene][t1_taxon]))
								genebinids_t1 = ','.join(temp)
								contig_taxa_connection = t1contigpres
								#contigtaxa1_t1 = 'contigtaxa_'+','.join(t1contigtaxa)
								contigfulltaxa1_t1= 'contigtaxa_'+','.join(t1contigfulltaxa)
								ho = has_higher_order(fulltaxa2,t1clusterfulltaxa)
								t1ll = lowest_level(t1clusterfulltaxa)
								counter=str(count)
								count+=1
								temp = []
								for t1_taxon in t1clustertaxa:
									temp.append(str(num_bincontigs[t1][level][gene][t1_taxon]))
								bincontigs_num_t1 = ','.join(temp)
								if bincontigs_num_t1 == ',':
									bincontigs_num_t1 = 'NA'
								bincontigs_num_t2 = str(num_bincontigs[t2][level][gene][taxa2])
								outfile2.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}\t{30}\t{31}\t{32}\n'.format(counter,connection_count,patient1,t1,t2,level,taxa1,taxa1_fulltaxa,ho,taxa1_t1_abund,taxa1_t2_abund,taxa2,fulltaxa2,taxa2_t1_abund,taxa2_t2_abund, gene,genename,genecontigs_t2,t1rpkm,t2rpkm,trans_taxa2_t1,trans_taxa2_t2,num_nont1,other_tps,contig_taxa_connection,contigfulltaxa1_t1,t1ll,trans_taxa1_t1,trans_taxa1_t2,bincontigs_num_t1,bincontigs_num_t2,genebinids_t1,genebinids_t2))
