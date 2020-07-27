import glob
import operator
import re
import sys

folder = sys.argv[1]
#doesn't handle shared most commons
def most_common(lst):
    return max(set(lst), key=lst.count)

pfamiddesc = {}
pfaminfo = open('/workdir/users/agk85/CDC2/mobile/metagenomes/pfam_acc_desc.txt')
for line in pfaminfo:
	pfamid = line.split('\t')[0].split('ACC   ')[1]
	desc = line.split('\t')[1].split('DESC  ')[1].strip()
	pfamiddesc[pfamid] = desc


pfamdict = {}
pfamiddict = {}
pfamdescdict = {}
pfampaths = glob.glob('/workdir/users/agk85/' + folder + '/annotation/metagenomes/*_pfam.txt')
for pfamfile in pfampaths:
	with open(pfamfile) as p:
		for line in p:
			pfamdict[line.split('\t')[0]] = line.split('\t')[2]
			pfamid = line.split('\t')[3]
			pfamiddict[line.split('\t')[0]] = pfamid
			desc = pfamiddesc[pfamid]
			pfamdescdict[line.split('\t')[0]] = desc

#mobile pfam genes
plasmid_pfams = []
with open('/workdir/users/agk85/CDC/annotation/Plasmid_pfamids.txt','r') as pfam_file:
	for line in pfam_file:
		plasmid_pfams.append(line.strip())

phage_pfams = []
with open('/workdir/users/agk85/CDC/annotation/Phage_pfamids.txt','r') as pfam_file:
	for line in pfam_file:
		phage_pfams.append(line.strip())

transposon_pfams = []
with open('/workdir/users/agk85/CDC/annotation/Transposon_pfamids.txt','r') as pfam_file:
	for line in pfam_file:
		transposon_pfams.append(line.strip())

other_pfams = []
with open('/workdir/users/agk85/CDC/annotation/Other_pfamids.txt','r') as pfam_file:
	for line in pfam_file:
		other_pfams.append(line.strip())


mobilepfamsource = {}
pfampaths = glob.glob('/workdir/users/agk85/' + folder + '/annotation/metagenomes/*_mobilegenes.txt')
for pfamfile in pfampaths:
	with open(pfamfile) as p:
		for line in p:
			tpn = 0
			pmd = 0
			phg = 0
			otr = 0
			gene = line.split('\t')[0]
			pfamid = line.split('\t')[7]
			if pfamid in transposon_pfams:
				tpn = 1
			if pfamid in plasmid_pfams:
				pmd = 1
			if pfamid in phage_pfams:
				phg = 1
			if pfamid in other_pfams:
				otr = 1
			mobilepfamsource[gene] = (pmd, phg, tpn, otr)

#gene sources
#phage|plasmid|transposon

mge_source = {}

#I grepped the '>' grep '>' mobile.fna > gene_ids.txt
with open('/workdir/users/agk85/' + folder + '/mobile/metagenomes/gene_ids.txt') as mobile:
	for line in mobile:
		gene = line.split(' ')[0].split('>')[1]
		plasmid = int(line.split(' ')[1].strip().split('|')[0])
		phage = int(line.split(' ')[1].strip().split('|')[1])
		transposon = int(line.split(' ')[1].strip().split('|')[2])
		ice = int(line.split(' ')[1].strip().split('|')[3])
		ime = int(line.split(' ')[1].strip().split('|')[4])
		arg = int(line.split(' ')[1].strip().split('|')[5])
		mobilepfam = int(line.split(' ')[1].strip().split('|')[6])
		other = 0
		if mobilepfam == 1:
			pmd, phg, tpn, otr = mobilepfamsource[gene]
			plasmid = plasmid + pmd
			phage = phage + phg
			transposon = transposon + tpn
			other = other + otr
		try:
			print(mge_source[gene])
		except KeyError:
			mge_source[gene] = (plasmid, phage, transposon, ice,ime,other,arg)

def binary(s):
	if s>0:
		b = 1
	else:
		b = 0
	return b

def mge_type(s1,s2,s3,s4,s5,s6):
	#we don't care about args
	b1 = binary(s1)
	b2 = binary(s2)
	b3 = binary(s3)
	b4 = binary(s4)
	b5 = binary(s5)
	b6 = binary(s6)
	mgetype = ''
	mgecombotype = ''
	splasmid = binary(s1 + s4) # add plasmids + ices
	stransposon =binary(s3 + s5) #add transposons + imes
	if (b1+b2+b3+b4+b5+b6)>1:
		mgetype = 'multiple'
	elif b1>0:
		mgetype = 'plasmid'
	elif b2>0:
		mgetype = 'phage'
	elif b3>0:
		mgetype = 'transposon'
	elif b4>0:
		mgetype = 'ice'
	elif b5>0:
		mgetype = 'ime'
	elif b6>0:
		mgetype = 'other'
	if (splasmid+b2+stransposon+s6)>1:
		mgecombotype = 'multiple'
	elif splasmid>0:
		mgecombotype = 'plasmid'
	elif b2>0:
		mgecombotype = 'phage'
	elif stransposon>0:
		mgecombotype = 'transposon'
	elif s6>0:
		mgecombotype = 'other'
	return mgetype,mgecombotype

#search a list of descriptions and return whether or not the search terms are in it for each of the mge types
transposon_terms = ['ranspos', 'nsertion element', 'is element', 'IS element','IS[0-9]']
phage_terms = ['hage', 'ail protein', 'tegument', 'apsid', 'ail fibre', 'ail assembly', 'ail sheath', 'ail tube']
plasmid_terms = ['elaxase', 'onjug', 'Trb', 'ype IV', 'ype iv', 'mob','Mob','Tra[A-Z]', 't4ss', 'T4SS', 'esolvase','plasmid','Vir[A-Z][0-9]']
other_terms =['ntegrase']


mgetypes = [plasmid_terms, phage_terms,transposon_terms, other_terms]
mgenames = ['plasmid','phage','transposon','other']

def search_descriptions(descriptions):
	"""goal is to get a code 0,1,1,0 telling if machinery terms in the descriptions"""
	code = ['0','0','0','0']
	for i in range(len(mgetypes)):
		mgetype = mgetypes[i]
		for term in mgetype:
			if re.search(term, descriptions):
				code[i] = '1'
	return code	

inhandle = '/workdir/users/agk85/' + folder + '/mobile/metagenomes/mge_99_nr.fna.clstr.tbl'
outhandle = '/workdir/users/agk85/' + folder + '/mobile/metagenomes/mge_99_nr.fna.clstr.tbl.annot'
deschandle = '/workdir/users/agk85/' + folder + '/mobile/metagenomes/mge_99_nr.fna.clstr.desc'
with open(outhandle,'w') as outfile:
	with open(deschandle,'w') as descfile:
		outfile.write('Cluster\tPfam_annotation\tMGE_type\tMGE_combo_type\tPlasmid\tPhage\tTransposon\tICE\tIME\tOther\tARG\tCluster2\tRep_gene\tGenes\tGene_count\tPfam_ids\n')
		descfile.write('Cluster\tPfam_descriptions\tPlasmid_mach\tPhage_mach\tTransposon_mach\tOther_mach\n')
		with open(inhandle) as infile:
			for line in infile:
				cluster = line.split('\t')[0]
				genenames = line.split('\t')[2].strip().split(',')
				representative = line.split('\t')[1]
				gene_annotations = []
				pfam_descriptions = []
				pfam_ids = []
				gene_sources = (0,0,0,0,0,0,0)
				for gene in genenames:
					try:
						gene_annot = pfamdict[gene]
						pfam_ids.append(pfamiddict[gene])
						pfam_descriptions.append(pfamdescdict[gene])
						gene_annotations.append(gene_annot)
					except KeyError:
						a = 1
					try:
						new_source = mge_source[gene]
						gene_sources = tuple(map(operator.add, gene_sources, new_source))
					except KeyError:
						a = 1
				if len(set(gene_annotations)) > 1:
					print(cluster)
					print(gene_annotations)
				try:
					rep_annot = pfamdict[representative]
					rep_source = mge_source[representative] 
				except KeyError:
					try:
						rep_annot = most_common(gene_annotations)
					except ValueError:
						rep_annot = 'NA'
				pfamids = ','.join(list(set(pfam_ids)))
				pfamdescriptions = ','.join(list(set(pfam_descriptions)))
				s1,s2,s3,s4,s5,s6,s7 = gene_sources
				mge_finaltype, mge_finalcombotype = mge_type(s1,s2,s3,s4,s5,s6)
				outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n'.format(cluster, rep_annot, mge_finaltype,mge_finalcombotype,s1, s2,s3,s4,s5,s6,s7,line.strip(), pfamids))
				descfile.write('{0}\t{1}\t{2}\n'.format(cluster, pfamdescriptions, '\t'.join(search_descriptions(pfamdescriptions))))



#afterwards you have to do this:
#sed -i -e "s/'/ prime/g" mge_95_nr.fna.clstr.desc  to make sure the primes are not disrupting R

