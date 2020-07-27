#contig_vs_bin.py
#USAGE:
#python contig_vs_bin.py -c ~/agk/press2/das/ProxiMeta-1/ProxiMeta-1_DASTool_scaffolds2bin.txt -o ~/agk/press2/bins/ProxiMeta-1_MluCI-1_das_5_contigs_v_bins.txt -l ~/agk/press2/hicpro/single/MluCI-1_output/hic_results/data/MluCI-1/MluCI-1_allValidPairs -t ~/agk/press2/combo_tables/metagenomes/ProxiMeta-1_master_scf_table.txt -m 5
#python contig_vs_bin.py -b /workdir/users/agk85/simulated/das/sim-1/sim-1_DASTool_scaffolds2bin.txt -o /workdir/users/agk85/simulated/bins/sim-1_das_5_contigs_v_bins.txt -l /workdir/users/agk85/simulated/hicpro/output/sim-1_output/hic_results/data/sim-1/sim-1_trans_primary_0_ncol_withexcise_noeuks_normalize.txt -t /workdir/users/agk85/simulated/combo_tables/metagenomes/sim-1_master_scf_table.txt -m 5
#-c sim-1_mge_count.qa


#WRK=/workdir/users/agk85/CDC2
#
# python contig_vs_bin.py \
# -b $WRK/das/B357-3/B357-3_DASTool_scaffolds2bin.txt \
# -o $WRK/bins/B357-3_das_2_contigs_v_bins.txt \
# -l $WRK/hicpro/output/B357-3_output/hic_results/data/B357-3/B357-3_trans_primary_0_ncol_withexcise_noeuks_normalize_2.txt \
# -t $WRK/combo_tables/metagenomes/B357-3_master_scf_table.txt \
# -m 2 -c $WRK/das/B357-3/checkm_lineage/B357-3.stats \
# -a $WRK/args/args_99_nr.fna.clstr.tbl \
# -k $WRK/das/B357-3/kraken/B357-3_all_kraken.txt \
# -kc $WRK/kraken/B357-3/B357-3.kraken.taxonomy.txt

import argparse
from argparse import RawDescriptionHelpFormatter
from copy import deepcopy

#def getOptions():
description="""Program creates a contig vs. bin table summing all contacts between the contig and the bin of interest, not counting contacts to itself (Only trans)"""

parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-b','--bins', dest="binhandle", action='store', required=True,  help='bin-contig file [REQUIRED]', metavar="INFILE")
parser.add_argument('-l','--hic', dest="hichandle", action='store', required=True, help='HI-C file [REQUIRED]', metavar="HIC_FILE")
parser.add_argument('-t','--table', dest="combotable", action='store', required=True, help='ComboTable file [REQUIRED]', metavar="COMBO_FILE")
parser.add_argument('-m','--min', dest="min", action='store', type=int, required=True, help='Min contacts [REQUIRED]', metavar="MIN_PROP")
parser.add_argument('-o','--out', dest="outhandle", action='store', required=True, help='Outfile [REQUIRED]', metavar="OUTFILE")
parser.add_argument('-c','--checkm', dest="checkm", action='store', required=True, help='checkm_info [REQUIRED]', metavar="CHECKM")
parser.add_argument('-ar','--argclust', dest="argclust", action='store', required=True, help='ARG clusters [REQUIRED]', metavar="ARGCLUST")
parser.add_argument('-mge','--mgeclust', dest="mgeclust", action='store', required=True, help='MGE clusters [REQUIRED]', metavar="MGECLUST")
parser.add_argument('-k','--kraken', dest="kraken_bins", action='store', required=False, help='Kraken cluster file', metavar="KRAKEN")
parser.add_argument('-kc','--krakenc', dest="kraken_contigs", action='store', required=True, help='Kraken contigs file ', metavar="KRAKEN_CONTIG")
#parser.add_argument('-g','--gtdb', dest="gtdb_clusters", action='store', required=True, help='GTDB cluster taxonomies [REQUIRED]', metavar="GTDB")


args = parser.parse_args()

sample = args.kraken_contigs.split('/')[-1].split('.kraken.taxonomy.txt')[0]
cluster_arg_dict = {}
cluster_mge_dict = {}
arg_cluster_dict = {} #so still deciding if I want the clusters OR if I want the args themselves---let's go with both!
mge_cluster_dict = {}
arg_gene_dict = {}
mge_gene_dict = {}



scf_argpresence={}
scf_argclusters={}
scf_args={}
#initialize the argdict
with open(args.argclust) as argfile:
	for line in argfile:
		cluster = line.split('\t')[0]
		arglist = line.split('\t')[2].split(',')
		cluster_arg_dict[cluster] = arglist
		for arg in arglist:
			scf = '_'.join(arg.split('_')[0:-1])
			arg_cluster_dict[arg] = cluster
			scf_argpresence[scf] = 1
			try:
				scf_argclusters[scf].append(cluster)
				scf_args[scf].append(arg)
			except KeyError:
				scf_argclusters[scf] = [cluster]
				scf_args[scf] = [arg]


scf_mgepresence={}
scf_mgeclusters={}
scf_mges={}
#initialize the argdict
with open(args.mgeclust) as mgefile:
	for line in mgefile:
		cluster = line.split('\t')[0]
		mgelist = line.split('\t')[2].split(',')
		cluster_mge_dict[cluster] = mgelist
		for mge in mgelist:
			scf = '_'.join(mge.split('_')[0:-1])
			mge_cluster_dict[mge] = cluster
			scf_mgepresence[scf] = 1
			try:
				scf_mgeclusters[scf].append(cluster)
				scf_mges[scf].append(mge)
			except KeyError:
				scf_mgeclusters[scf] = [cluster]
				scf_mges[scf] = [mge]

#initialize the hic individual dictionaries
hic_count_individual = {}
hic_norml_individual = {}
hic_normrf_individual = {}
contigs = []

length_dict = {}
#get all the contigs
with open(args.combotable) as infile:
	header = infile.readline()
	for line in infile:
		contig = line.split('\t')[0]
		length = line.split('\t')[1].split('_')[1]
		length_dict[contig] = length
		contigs.append(contig)
		hic_count_individual[contig] = 0
		hic_norml_individual[contig] = 0
		hic_normrf_individual[contig] = 0


#checkm quality
qualitydict = {}
cluster_length = {}
with open(args.checkm) as checkm:
	for line in checkm:
		if 'network' in args.outhandle:
			bin = line.split('\t')[1].split('.contigs')[0]
			completion = float(line.split('\t')[3])
			contamination = float(line.split('\t')[4])
			length = line.strip().split('\t')[6]
		else:
			bin = line.split('\t')[0].split('.contigs')[0]
			completion = float(line.split('\t')[2])
			contamination = float(line.split('\t')[3])
			length = line.strip().split('\t')[5]
			print(completion, contamination)
		quality = 'BAD'
		if (completion > 90 and contamination < 5):
			quality = 'HQ'
		elif (completion >= 50 and contamination < 10):
			quality = 'MQ'
		elif (completion < 50 and contamination < 10):
			quality = 'LQ'	
		qualitydict[bin] = quality
		cluster_length[bin] = length


#start with full list and remove iteratively as you find them in the bins
#do not include bad bins
noncluster_contigs = deepcopy(contigs)
####BIN STUFF######
#create a dictionary of the contigs and bins
binlist = []
bincontigsdict = {} #put in bin and get contigs
contigbindict = {} #put in contig and get bin (meaning that some will not have a bin)
with open(args.binhandle) as binfile:
	for line in binfile:
		bin = line.strip().split('\t')[1]
		if qualitydict[bin] != 'BAD':
			binlist.append(bin) #this should get all of them
			contig = line.split('\t')[0]
			try:
				bincontigsdict[bin].append(contig)		
			except KeyError:
				bincontigsdict[bin] = [contig]
			contigbindict[contig] = bin
			#get rid of this contig (because it has a bin)
			try:
				noncluster_contigs.remove(contig)
			except:
				a =1 #it's a kraken euk!

bins = list(set(binlist))



print(qualitydict)
#kraken bin taxonomies
#B357-3.054_sub.contigs.fa.report.txt.besttaxid  75.0    2       superkingdom    1       k__Bacteria; p__; c__; o__; f__; g__; s__;
	
kraken_bin_taxonomy = {}
kraken_bin_contig_taxonomy = {}
if 'das' in args.kraken_bins or 'bin3c' in args.kraken_bins:
	with open(args.kraken_bins) as krakenfile:
		for line in krakenfile:
			bin = line.split('.contigs.fa.report.txt.besttaxid')[0]
			try:
				contiglist = bincontigsdict[bin]
				taxonomy = line.strip().split('\t')[7]
				kraken_bin_taxonomy[bin] = taxonomy
				for contig in contiglist:
					kraken_bin_contig_taxonomy[contig] = taxonomy
			except KeyError:
				a = 1

if 'network' in args.kraken_bins:
	with open(args.kraken_bins) as krakenfile:
		header = krakenfile.readline()
		for line in krakenfile:
			bin = line.split('\t')[1]
			taxonomy = line.strip().split('\t')[3]
			kraken_bin_taxonomy[bin] = taxonomy
			contiglist = bincontigsdict[bin]
			for contig in contiglist:
				kraken_bin_contig_taxonomy[contig] = taxonomy

kraken_contig_taxonomy = {}
with open(args.kraken_contigs) as kraken:		
	kraken_dict = {}
	markerlist = ['d','p','c','o','f','g','s']
	for line in kraken:
			base_marker = {'d':'','p':'','c':'','o':'','f':'','g':'','s':''}
			scfid = line.split('\t')[0]
			taxonomy = line.strip().split('\t')[1]
			taxonomy_fields = taxonomy.split('|')
			for field in taxonomy_fields:
					if field != 'root':
							marker = field.split('__')[0]
							taxa = field.split('__')[1]
							if marker in markerlist:
									base_marker[marker] = taxa
			parsedtaxonomy = ''
			for marker in markerlist:
					#standard marker is 'k' for kingdom instead of domain, so switching that
					if marker == 'd':
							parsedtaxonomy = parsedtaxonomy + 'k__' + base_marker[marker]
					else:
							parsedtaxonomy = parsedtaxonomy + marker + '__' + base_marker[marker]
					#my taxonomy format stupidly has semicolon then spaces except after the species
					if marker != 's':
							parsedtaxonomy = parsedtaxonomy + '; '
					else:
							parsedtaxonomy = parsedtaxonomy + ';'
			kraken_contig_taxonomy[scfid] = parsedtaxonomy


#need to redo once rerunning network and checkm
# gtdb_bin_taxonomy = {}
# with open(args.gtdb_bins) as gtdbfile:
# 	for line in gtdbfile:
# 		bin = line.split('\t')[0]
# 		taxonomy = line.split('\t')[1]
# 		#convert to my format
# 		taxa1 = taxonomy.replace('d__', 'k__')
# 		taxa2 = taxa1.replace(' ', '_')
# 		taxa3 = taxa2.replace(';', '; ')
# 		taxa4 = taxa3 + ';'
# 		gtdb_bin_taxonomy[bin] = taxa4 #although there is something to be said about keeping it in their format!
# 

#get hic trans contig information from normalized info
#use the normalized file with counts column
#SRR6131122.10558730     ProxiMeta-1_1|phage|210862|217963_98338 110     -       ProxiMeta-1_scaffold_1  210442  +       382     HIC_ProxiMeta-1_1|phage|210862|217963_98338_5   HIC_ProxiMeta-1_scaffold_1_1823 42      42
hic_count_dict = {}
hic_norml_dict = {}
hic_normrf_dict = {}
with open(args.hichandle) as hicfile:
	for line in hicfile:
		#Get the contigs linking
		contig1 = line.split('\t')[0]
		contig2 = line.split('\t')[1]
		#make sure that it is a trans read!
		if contig1 != contig2:
			norml = float(line.split('\t')[2])
			normrf = float(line.split('\t')[3])
			count = float(line.split('\t')[4])
			try:
				hic_count_dict[contig1+contig2]+= count
				hic_count_dict[contig2+contig1]+= count
				hic_norml_dict[contig1+contig2]+= norml
				hic_norml_dict[contig2+contig1]+= norml
				hic_normrf_dict[contig1+contig2]+= normrf
				hic_normrf_dict[contig2+contig1]+= normrf
			except KeyError:
				hic_count_dict[contig1+contig2]= count
				hic_count_dict[contig2+contig1]= count
				hic_norml_dict[contig1+contig2]= norml
				hic_norml_dict[contig2+contig1]= norml
				hic_normrf_dict[contig1+contig2]= normrf
				hic_normrf_dict[contig2+contig1]= normrf
			#update all of the individuals with the counts
			hic_count_individual[contig1]+=count
			hic_norml_individual[contig1]+=norml
			hic_normrf_individual[contig1]+=normrf
			hic_count_individual[contig2]+=count
			hic_norml_individual[contig2]+=norml
			hic_normrf_individual[contig2]+=normrf




header = 'count\tsample\tcontig\tbin\tbin_length\tquality\tassociation_type\ttotal_count\ttrans_count\tnorm_l\tnorm_rf\tcontig_taxonomy\tkraken_taxonomy\tgtdb_taxonomy\targ_presence\targ_clusters\targ_genes\tmge_presence\tmge_clusters\tmge_genes\n'
#count so you have unique rownames
#association_type = cluster_resident, contig_resident, cluster_hic, contig_hic | cluster_resident overrides cluster_hic hic just means they are linked only by hic
#counts of trans, norm_l, and norm_rf are summed for the members in the cluster
#get the kraken and gtdb taxonomies ---# assoc. = # any entity, # any cluster, # unique taxa, # unique cluster taxa
#do i include the contig taxa or no, maybe the thing is ot do this 

#so either you are in it and it prints no matter what
#or you aren't in it and it prints the contigs 
#I don't want to add fake counts so we will ahve to make sure things are prepared beforehand


count = 0

#create contig table with hic data (long format)
with open(args.outhandle, 'w') as outfile:
	outfile.write(header)
	for contig1 in contigs:
		#1#########go through all of the clusters/non-clustered contigs and see if there in the bins
		for bin in bins:
			quality = 'NA' #until you update this
			pass_threshold = 0
			contigcount = 0 #count summed for the particular cluster/noncluster_contig
			contignorml = 0
			contignormrf = 0
			quality = qualitydict[bin]
			binlength = cluster_length[bin]
			for contig2 in bincontigsdict[bin]:
				if contig1 != contig2:
					try:
						contigcount += hic_count_dict[contig1+contig2]
						contignorml += hic_norml_dict[contig1+contig2]
						contignormrf += hic_normrf_dict[contig1+contig2]
					except KeyError:
						pass
			if contig1 in bincontigsdict[bin]:
				association_type = 'cluster_resident'
				pass_threshold = 1
			else:
				association_type = 'cluster_hic'
			#if float(contigcount/contig_total_count)>args.minprop: #threshold on percent of the cluster
			if contigcount>args.min:
				pass_threshold = 1
			if pass_threshold == 1:
				try:
					kraken_taxonomy = kraken_bin_taxonomy[bin] 
				except KeyError:
					kraken_taxonomy = 'NA'
				contig_total_count = hic_count_individual[contig1]
				#see if there are any args
				try:
					arg_present = scf_argpresence[contig1]
					arg_clusters = ','.join(scf_argclusters[contig1])
					arg_genes = ','.join(scf_args[contig1])
				except KeyError:
					arg_present = '0'
					arg_clusters = 'NA'
					arg_genes = 'NA'
				#see if there are any mges'
				try:
					mge_present = scf_mgepresence[contig1]
					mge_clusters = ','.join(scf_mgeclusters[contig1])
					mge_genes = ','.join(scf_mges[contig1])
				except KeyError:
					mge_present = '0'
					mge_clusters = 'NA'
					mge_genes = 'NA'
				#I'll run this later
				gtdb_taxonomy = 'NA'			
				try:
					#it should be the bin first, and if it doesn't have a taxonomy, then use the contig taxonomy
					contig_kraken_taxonomy = kraken_bin_contig_taxonomy[contig1] 
				except KeyError:
					try:
						contig_kraken_taxonomy = kraken_contig_taxonomy[contig1]
					except KeyError:
						contig_kraken_taxonomy = 'NA'
				outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\n'.format(str(count), sample, contig1, bin, binlength, quality, association_type,str(contig_total_count), str(contigcount), str(contignorml), str(contignormrf), contig_kraken_taxonomy, kraken_taxonomy, gtdb_taxonomy, arg_present, arg_clusters, arg_genes, mge_present, mge_clusters, mge_genes))
				count +=1
		#2############go through the rest of the contigs 
		for contig2 in noncluster_contigs:
			quality = 'NA'
			pass_threshold = 0
			contigcount = 0
			contignorml = 0
			contignormrf = 0
			try:
				contigcount += hic_count_dict[contig1+contig2]
				contignorml += hic_norml_dict[contig1+contig2]
				contignormrf += hic_normrf_dict[contig1+contig2]
			except KeyError:
				pass
			if contig1 == contig2:
				pass_threshold = 1
				association_type = 'contig_resident'
			else:
				association_type = 'contig_hic'
			if contigcount>1:
				#it's possible the contig is itself, and we still want to print it
				pass_threshold = 1
			if pass_threshold == 1:
				try:
					kraken_taxonomy = kraken_contig_taxonomy[contig2]
				except KeyError:
					kraken_taxonomy = 'NA'
				contig_total_count = hic_count_individual[contig1]
				try:
					arg_present = scf_argpresence[contig1]
					arg_clusters = ','.join(scf_argclusters[contig1])
					arg_genes = ','.join(scf_args[contig1])
				except KeyError:
					arg_present = '0'
					arg_clusters = 'NA'
					arg_genes = 'NA'
				#see if there are any mges'
				try:
					mge_present = scf_mgepresence[contig1]
					mge_clusters = ','.join(scf_mgeclusters[contig1])
					mge_genes = ','.join(scf_mges[contig1])
				except KeyError:
					mge_present = '0'
					mge_clusters = 'NA'
					mge_genes = 'NA'
				#I'll run this later
				gtdb_taxonomy = 'NA'	
				binlength = length_dict[contig2]
				try:
					contig_kraken_taxonomy = kraken_bin_contig_taxonomy[contig1] 
				except KeyError:
					try:
						contig_kraken_taxonomy = kraken_contig_taxonomy[contig1]
					except KeyError:
						contig_kraken_taxonomy = 'NA'
				outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\n'.format(str(count), sample,contig1, contig2, binlength, quality, association_type,str(contig_total_count), str(contigcount), str(contignorml), str(contignormrf), contig_kraken_taxonomy, kraken_taxonomy, gtdb_taxonomy, arg_present, arg_clusters, arg_genes, mge_present, mge_clusters, mge_genes))
				count +=1

#ok think about what you want to plot
#for every contig, how many OTHER taxonomies is it associated with, that aren't it
#so...kraken_taxonomy has to be from bin or contig2
#but do we need to know contig1's affiliation too...i think so




