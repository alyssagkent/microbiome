#phage comparison to known data

#1 There was no known hit., we assign 1 taxa. (putative host-specific phage)
#2 There was no known hit, we assign multiple taxa (putative broad-range phage)

#3 There was 1 hit, we get the same (confirms host-specific phage)
#4 There was 1 hit, we get 1 different hit (putative broad range, unknown before, or crappy hit)
#5 There was 1 hit, we get more hits, including the BLAST hit (putative broad range, unknown before, or crappy hit)
#5_2 There was 1 hit, we get more hits, not including the BLAST hit (putative broad range, unknown before, or crappy hit)

#6 There was multiple hits, we get some overlap (confirms at least a bt=it of what we know)
#7 There was multiple hits, we get 1 or more, none of which overlap (expand our knowledge of organisms it can invade).

#8 There were multiple hits and one hictaxa that was in it
#9 There were multiple hits and one hictaxa that was not in it

from taxid_to_taxonomy import taxid_dict
import sys

#get taxid to taxonomy per rankedlineage.dmp and node.dmp 
taxid_to_taxonomy = taxid_dict()

# fi there are bad taxids then you have to find them online and put them in your format
# taxid_to_taxonomy['2024339']='k__Viruses; p__; c__; o__Caudovirales; f__Siphoviridae; g__C2virus; s__Lactococcus_phage_62403;'

minreads = sys.argv[1]

bad_taxids = []
inhandle = '/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.out'
gene_taxonomies = {}
with open(inhandle) as infile:
	for line in infile:
		gene = line.split('\t')[0]
		taxids = line.strip().split('\t')[8].split(';')
		for taxid in taxids:
			try:
				taxonomy = taxid_to_taxonomy[taxid]
				try:
					gene_taxonomies[gene].append(taxonomy)
				except KeyError:
					gene_taxonomies[gene]=[taxonomy]
			except KeyError:
				bad_taxids.append(taxid)

print(set(bad_taxids))

cluster_hictaxa = {}
cluster_ncbitaxa = {}
cluster_basetaxa = {}
genedict = {}
noblastresult = []
comparehandle = '/workdir/users/agk85/CDC2/bins/patric_figures/mge_' + minreads + '_base_hic_taxonomies.txt'
#Cluster Gene    Source  Number_base_taxa        Number_hic_taxa Number_base_species     Number_hic_species      Base_taxonomies Hic_taxonomies  Base_species    Hic_species     Pfams
with open(comparehandle) as cfile:
	header = cfile.readline()
	for line in cfile:
		cluster = line.split('\t')[0]
		gene = line.split('\t')[1]
		try:
			genedict[cluster].append(gene)
		except KeyError:
			genedict[cluster] = [gene]
		hic_taxonomies = line.split('\t')[8].split(',')
		base_taxonomies = line.split('\t')[7].split(',')
		hicsupport_taxonomies = line.split('\t')[10].strip().split(',')
		hic_taxonomies = list(set(hic_taxonomies + hicsupport_taxonomies))
		if len(hic_taxonomies)>1:
			try:
				hic_taxonomies.remove('')
			except ValueError:
				a=1
		if hic_taxonomies != ['']:
			try:
				cluster_hictaxa[cluster] = cluster_hictaxa[cluster] + hic_taxonomies
			except KeyError:
				cluster_hictaxa[cluster] = hic_taxonomies
			if base_taxonomies != ['']:
				try:
					cluster_basetaxa[cluster] = cluster_basetaxa[cluster] + base_taxonomies
				except KeyError:
					cluster_basetaxa[cluster] = base_taxonomies
			try:
				cluster_ncbitaxa[cluster] = cluster_ncbitaxa[cluster] + gene_taxonomies[gene]
			except KeyError:
				try:
					cluster_ncbitaxa[cluster] = gene_taxonomies[gene]
				except:
					#blast result did not return anything
					noblastresult.append(gene)


len(noblastresult) #2579  | with hic taxa | 1221
len(cluster_hictaxa) #2587 | with hic taxa | 1079
len(cluster_ncbitaxa) #1479 | with hic taxa | 647
len(cluster_basetaxa) 

def highertaxonomies(taxalist, delim):
	newtaxalist = []
	for taxonomy in taxalist:
		leveltaxa = taxonomy.split(delim)[1].split(';')[0]
		if (leveltaxa != '' and 'k__Viruses' not in taxonomy):
			newtaxonomy = taxonomy.split(delim)[0]+delim + leveltaxa + ';'
			newtaxalist.append(newtaxonomy)
	newtaxasetlist = list(set(newtaxalist))
	return newtaxasetlist


def notvirustaxonomies(taxalist):
	newtaxalist = []
	for taxonomy in taxalist:
		if ('k__Virus' not in taxonomy):
			newtaxalist.append(taxonomy)
	return newtaxalist

#iterate over the hic clusters
#what number of them 
delims = ['s__','g__','f__','o__','c__','p__','k__']
delimnames = ['species','genus','family','order','class','phylum','kingdom']
#delims = ['s__']
outhandle = '/workdir/users/agk85/CDC2/bins/patric_figures/mge_' + minreads + '_support_ncbi_distribution.txt'
with open(outhandle,'w') as outfile:
	outfile.write('delimname\tnum1\tnum2\tnum3\tnum4\tnum5\tnum5_2\tnum6\tnum7\tnum8\tnum8_2\tnum9\tunannotated_level\ttotal\tnohit\toverlap\tno_overlap\texpand\n')
	for i in range(len(delims)):
		num1 = 0
		num2 = 0
		num3 = 0
		num4 = 0
		num5 = 0
		num5_2 = 0
		num6 = 0
		num7 = 0
		num8 = 0
		num8_2 = 0
		num9 = 0
		missing_hicname = 0
		missing_ncbiname = 0
		unannotated_level = 0
		idk = 0
		delim = delims[i]
		delimname = delimnames[i]
		for cluster in cluster_hictaxa.keys():
			hictaxa = highertaxonomies(cluster_hictaxa[cluster],delim)
			orighic = cluster_hictaxa[cluster]
			if cluster in cluster_ncbitaxa.keys():
				ncbitaxa = highertaxonomies(cluster_ncbitaxa[cluster],delim)
				origncbi = cluster_ncbitaxa[cluster]
				notvirusncbi = notvirustaxonomies(cluster_ncbitaxa[cluster])
				#if 'k__Virus' in ''.join(cluster_ncbitaxa[cluster]):
					#print cluster
					#print cluster_ncbitaxa[cluster]
			else:
				ncbitaxa = []
				origncbi =  []
				notvirusncbi = []
			if (len(hictaxa) == 0):
				missing_hicname = missing_hicname + 1
				unannotated_level = unannotated_level + 1
			elif (len(hictaxa)==1):
				if (len(ncbitaxa) == 0) and (len(origncbi)-len(notvirusncbi) ==0) and (len(origncbi)>0):
					unannotated_level = unannotated_level + 1
					missing_ncbiname = missing_ncbiname + 1
					#what about if empty and virus or empty and not just taxa unannotated:
				elif ((len(ncbitaxa) == 0) and (len(origncbi)-len(notvirusncbi) >0) and (len(origncbi)>0)) or ((len(ncbitaxa) == 0) and (len(origncbi)==0)):
					num1 = num1 + 1
				elif (len(ncbitaxa) == 1):
					if (ncbitaxa[0] == hictaxa[0]):
						#1 to 1 contained
						num3 = num3+1
					else:
						#1 to 1 not overlapped
						num4 = num4 + 1
				elif (len(ncbitaxa) > 1):
					if hictaxa[0] in ncbitaxa:
						num6 = num6 + 1
					else:
						num7 = num7 + 1
				else:
					idk = idk + 1
					print(cluster)
					print(origncbi)
					print(orighic)	
			elif (len(hictaxa)>1):
				if (len(ncbitaxa) == 0) and (len(origncbi)-len(notvirusncbi) ==0) and (len(origncbi)>0):
					unannotated_level = unannotated_level + 1
					missing_ncbiname = missing_ncbiname + 1
				#what about if empty and virus or empty and not just taxa unannotated:
				elif ((len(ncbitaxa) == 0) and (len(origncbi)-len(notvirusncbi) >0) and (len(origncbi)>0)) or ((len(ncbitaxa) == 0) and (len(origncbi)==0)):
					num2 = num2 + 1
				elif (len(ncbitaxa) == 1):
					#the single ncbitaxa is in the group of many hictaxa
					if ncbitaxa[0] in hictaxa:
						num5= num5 + 1
					#the single ncbitaxa is not in the many hictaxa
					else:
						num5_2=num5_2 + 1
				elif (len(ncbitaxa) > 1):
					flag = 0
					for taxa in ncbitaxa:
						if taxa in hictaxa:
							flag = 1
					if flag == 1:
						#do they completely overlap?
						if (len(set(hictaxa)-set(ncbitaxa))==0 ):	
							num8 = num8 + 1
						else:
							num8_2 = num8_2 + 1
					if flag == 0:
						num9 = num9 + 1
				else:
					idk = idk + 1
					print(cluster)
					print(origncbi)
					print(orighic)			
		nohit = num1 + num2
		overlap = num3 + num5 + num6 + num8_2 + num8
		nooverlap = num4 + num5_2 + num7 + num9
		expand = num1 + num2 + num4 + num5 + num5_2 + num8_2 + num7 + num9
		total = num1 + num2 + num3 + num4 + num5 + num5_2 + num6 +  num7 + num8 +num8_2 + num9 + unannotated_level
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\n'.format(delimname, num1, num2, num3, num4, num5, num5_2, num6, num7, num8, num8_2, num9, unannotated_level, total, nohit, overlap, nooverlap,expand))						


						
# print delim
# num1 #no ncbi and one hic
# num2 #no ncbi and multiple hic
# num3 #one hit one base
# num4 #one hit one different hic
# num5 # one hit and multiple hic including hit
# num5_2 #one hit and multiple hic not including hit
# num6 #multiple hits and all taxa overlap
# num6_2 #multiple hits and at least one overlaps
# num7 #multiple hits and none overlap
# num8 #multiple hits and multiple hic overlap
# num9 #multiple hits and multiple hic no overlap
# num1 + num2 + num3 + num4 + num5 + num5_2 + num6 + num6_2 + num7 + num8 + num9 + unannotated_level

