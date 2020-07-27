import glob
from taxid_to_taxonomy import taxid_dict

topargs = []
refhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/topargs/all_clusters_topargs_mergednames.txt'
with open(refhandle) as reffile:
	for line in reffile:
		cluster = line.split('\t')[0]
		topargs.append(cluster)


#get taxid to taxonomy per rankedlineage.dmp and node.dmp 
taxid_to_taxonomy = taxid_dict()

taxid_to_taxonomy['2024339']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_62403;'
taxid_to_taxonomy['2024338']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_62402;'
taxid_to_taxonomy['2024335']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_37203;'
taxid_to_taxonomy['2024334']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_05802;'
taxid_to_taxonomy['2024337']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_50504;'
taxid_to_taxonomy['2024336']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_50102;'
taxid_to_taxonomy['2282309']='k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__; s__Enterobacteriaceae_bacterium_w17;'
taxid_to_taxonomy['2282310']='k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__; s__Enterobacteriaceae_bacterium_w6;'
taxid_to_taxonomy['2268104']='k__other sequences; p__artificial sequences; c__vectors; o__; f__; g__; s__Phagemid_vector_pScaf-8064.3;'
taxid_to_taxonomy['2268096']='k__other sequences; p__artificial sequences; c__vectors; o__; f__; g__; s__Phagemid_vector_pScaf-5544.1;'
taxid_to_taxonomy['2268106']='k__other sequences; p__artificial sequences; c__vectors; o__; f__; g__; s__Phagemid_vector_pScaf-8064.5;'
taxid_to_taxonomy['2268101']='k__other sequences; p__artificial sequences; c__vectors; o__; f__; g__; s__Phagemid_vector_pScaf-7560.2;'
taxid_to_taxonomy['2268097']='k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillales Family XI. Incertae Sedis; g__Gemella; s__Phagemid_vector_pScaf-5544.2;'
taxid_to_taxonomy['2040624']='k__Viruses; p__; c__; o__Caudovirales; f__Myoviridae; g__Tevenvirinae; s__T4virusGemella_sp._ND_6198;'
taxid_to_taxonomy['2234085']='k__Viruses; p__; c__; o__Caudovirales; f__Myoviridae; g__Tevenvirinae; s__T4virusEscherichia_phage_vB_EcoM_JB75;'
taxid_to_taxonomy['2202142']='k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Chromobacteriaceae; g__Chromobacterium group; s__ChromobacteriumChromobacterium_sp._IIBBL_274-1;'
taxid_to_taxonomy['2234047']='k__Viruses; p__; c__; o__Caudovirales; f__Myoviridae; g__Tevenvirinae; s__Jd18virusKlebsiella_phage_Mineola;'
taxid_to_taxonomy['2268610']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__Tunavirinae; s__Kp36virusKlebsiella_phage_NJS1;'
taxid_to_taxonomy['2268098']='k__other sequences; p__artificial sequences; c__vectors; o__; f__; g__; s__Phagemid_vector_pScaf-5544.3;'
taxid_to_taxonomy['2268099']='k__other sequences; p__artificial sequences; c__vectors; o__; f__; g__; s__Phagemid_vector_pScaf-5544.4;'
taxid_to_taxonomy['2267236']='k__Viruses; p__; c__; o__Caudovirales; f__Myoviridae; g__Tevenvirinae; s__MoonvirusCitrobacter_phage_CF1_ERZ-2017;'
taxid_to_taxonomy['2024340']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_62606;'
taxid_to_taxonomy['2024341']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_74001;'


taxid_to_taxonomy['1739443']='k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Achromobacter; s__Achromobacter_sp._HMSC056C09;'
taxid_to_taxonomy['1581052']='k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Achromobacter; s__Achromobacter_sp._HMSC18C08;'
taxid_to_taxonomy['469608']='k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Klebsiella; s__Klebsiella_sp._1_1_55;'
taxid_to_taxonomy['469595']='k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Citrobacter; s__Citrobacter_sp._30_2;'
taxid_to_taxonomy['1739526']='k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Achromobacter; s__Achromobacter_sp._HMSC057D05;'
taxid_to_taxonomy['1581105']='k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Achromobacter; s__Achromobacter_sp._HMSC15D03;'
taxid_to_taxonomy['469598']='k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia_sp._3_2_53FAA;'

taxid_to_taxonomy['657310']='k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Enterococcus_sp._7L76;'
taxid_to_taxonomy['665937']='k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Lachnospiraceae; g__Anaerostipes; s__Anaerostipes_sp._3_2_56FAA;'
taxid_to_taxonomy['563191']='k__Bacteria; p__Firmicutes; c__Negativicutes; o__Acidaminococcales; f__Acidaminococcaceae; g__Acidaminococcus; s__Acidaminococcus_sp._D21;'
taxid_to_taxonomy['469594']='k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium; s__Bifidobacterium_sp._12_1_47BFAA;'
taxid_to_taxonomy['469596']='k__Bacteria; p__Firmicutes; c__Erysipelotrichia; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Coprobacillus; s__Coprobacillus_sp._29_1;'

taxid_to_taxonomy['469586']='k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_sp._1_1_6;'
taxid_to_taxonomy['665938']='k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_sp._2_1_56FAA;'
taxid_to_taxonomy['469587']='k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides_sp._2_1_16;'

taxid_to_taxonomy['702953'] ='k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Corynebacteriales; f__Corynebacteriaceae; g__Corynebacterium; s__Corynebacterium sp. NML00-0156;'
taxid_to_taxonomy['2211115'] ='k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Enterobacter; s__Enterobacter cloacae complex'
taxid_to_taxonomy['702969'] ='k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Corynebacteriales; f__Corynebacteriaceae; g__Corynebacterium; s__Corynebacterium sp. NML99-0020;'
taxid_to_taxonomy['702959'] ='k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Corynebacteriales; f__Corynebacteriaceae; g__Corynebacterium; s__Corynebacterium sp. NML93-0607;'
taxid_to_taxonomy['2250596'] ='k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus sp. 596553;'

pidthresh = 95
bad_taxids = []
inhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/patric_comparisons/toparg_prots_1e-100.out'
gene_taxonomies = {}
with open(inhandle) as infile:
	for line in infile:
		pid = float(line.split('\t')[3])
		if pid >= pidthresh:
			gene = line.split('\t')[1]
			taxid = line.strip().split('\t')[0]
			try:
				taxonomy = taxid_to_taxonomy[taxid]
				try:
					gene_taxonomies[gene].append(taxonomy)
				except KeyError:
					gene_taxonomies[gene]=[taxonomy]
			except KeyError:
				bad_taxids.append(taxid)

set(bad_taxids)

samples = []
samplist = glob.glob('/workdir/users/agk85/CDC/todo/*')
for samplong in samplist:
	samp = samplong.split('/')[6]
	samples.append(samp)

inhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/args_95_nr.fna.clstr'
protclusterdict = {}
clusterprotdict = {}
clusternum_map = {}
with open(inhandle) as infile:
	for line in infile:
		if line[0] == '>':
			cluster = line.strip().split(' ')[1]
		else:
			prot =  line.split('>')[1].split('...')[0]
			samp = prot.split('_')[0]
			protclusterdict[prot] = cluster
			if samp in samples:
				try:
					clusterprotdict[cluster].append(prot)
				except:
					clusterprotdict[cluster] = [prot]
			if '*' in line: 
				clusternum_map[cluster] = prot

cluster_ncbitaxa = {}
for cluster in clusterprotdict.keys():
	if cluster in topargs:
		genes = clusterprotdict[cluster]
		for gene in genes:
			sample = gene.split('_')[0]
			if sample in samples:
				try:
					cluster_ncbitaxa[cluster] = cluster_ncbitaxa[cluster] + gene_taxonomies[gene]
				except KeyError:
					try:
						cluster_ncbitaxa[cluster] = gene_taxonomies[gene]
					except:
						print ('does not have ncbi taxonomy')


outhandle = '/workdir/users/agk85/CDC/arg_v_org/metagenomes3/taxonomic_distributions/Patric_toparg_taxonomies.txt'
with open(outhandle,'w') as outfile:
	for cluster in cluster_ncbitaxa.keys():
		taxonomies = list(set(cluster_ncbitaxa[cluster]))
		for taxonomy in taxonomies:
			outfile.write('{0}\t{1}\n'.format(cluster, taxonomy))
