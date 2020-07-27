#USAGE: python get_taxonomy_overlap.py /workdir/users/agk85/press2/combo_tables/metagenomes/ProxiMeta-1_master_scf_table.txt /workdir/users/agk85/press2/combo_tables/metagenomes/ProxiMeta-1_taxonomy_overlap.txt ProxiMeta-1
#goal: get the shared taxa between gaemr kraken and besttaxa
import sys
inhandle = sys.argv[1]
outhandle = sys.argv[2]
sample = sys.argv[3]

keys = ['k','p','c','o','f','g','s']

kraken_dict ={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
gaemr_dict ={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
best_dict ={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
kg_shared_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
kb_shared_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
gb_shared_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
kgb_shared_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}

#annotation dictionaries---did it have any sort of annotation in that position?
kg_annotated_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
kb_annotated_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
gb_annotated_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
kgb_annotated_dict={'k':0, 'p':0, 'c':0, 'o':0, 'f':0, 'g':0, 's':0}
linecount = 0
with open(inhandle) as infile:
	header = infile.readline()
	for line in infile:
		linecount +=1
		kraken = line.split('\t')[35]
		gaemr = line.split('\t')[36].split('|')[0]
		best = line.strip().split('\t')[-1]
		#contig_length = float(line.split('length_')[1].split('\t')[0])
		for key in keys:
			if kraken != '.':
				krakentaxa = kraken.split(key + '__')[1].split(';')[0]
				krakentaxa = krakentaxa.replace(' ', '_')
			else: 
				krakentaxa = ''
			if (gaemr != '.' and gaemr != 'gaemr_.'):
				gaemrtaxa = gaemr.split(key + '__')[1].split(';')[0]
				gaemrtaxa = gaemrtaxa.replace(' ', '_')
			else:
				gaemrtaxa = ''
			if best != '.':
				besttaxa = best.split(key + '__')[1].split(';')[0]
				besttaxa.replace(' ', '_')
			else:
				besttaxa = ''
			if krakentaxa != '':
				kraken_dict[key] +=1
				#k + g annotated?)
				if gaemrtaxa != '':
					kg_annotated_dict[key] +=1
					#all three annotated?
					if besttaxa != '':
						kgb_annotated_dict[key] +=1
				#k + g the same
				if gaemrtaxa == krakentaxa:
					kg_shared_dict[key] +=1
					#all three shared?
					if besttaxa == krakentaxa:
						kgb_shared_dict[key] +=1
				#k + b annotated?
				if besttaxa!='':
					kb_annotated_dict[key] +=1
				#k + b the same?
				if besttaxa == krakentaxa:
					kb_shared_dict[key] +=1
			if gaemrtaxa != '':
				gaemr_dict[key] +=1
				#g + b annotated?
				if besttaxa!='':
					gb_annotated_dict[key] +=1
				#g + b the same?
				if besttaxa == gaemrtaxa:
					gb_shared_dict[key] +=1
			if besttaxa != '':
				best_dict[key] +=1


groups = ['kraken', 'gaemr','best', 'kg','kb','gb','kgb']
#link these groups to their respective dictionaries
groupdict = {'kraken':kraken_dict, 'gaemr':gaemr_dict,'best':best_dict, 'kg':kg_shared_dict,'kb':kb_shared_dict,'gb':gb_shared_dict,'kgb':kgb_shared_dict, 'kg_both_annotated':kg_annotated_dict,'kb_both_annotated':kb_annotated_dict,'gb_both_annotated':gb_annotated_dict,'kgb_both_annotated':kgb_annotated_dict}
#and to the dictionarys used for the denominators
annotateddict={'kraken':kraken_dict, 'gaemr':gaemr_dict,'best':best_dict,'kg':kg_annotated_dict,'kb':kb_annotated_dict,'gb':gb_annotated_dict,'kgb':kgb_annotated_dict}

header = 'Sample\tGroup\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tTotal_contigs\tKingdom_proportion\tPhylum_proportion\tClass_proportion\tOrder_proportion\tFamily_proportion\tGenus_proportion\tSpecies_proportion\n'
with open(outhandle,'w') as outfile:
        outfile.write(header)
        for group in groups:
		outfile.write(sample)
        	outfile.write('\t' + group)
        	for key in keys:
                	outfile.write('\t' + str(groupdict[group][key]))
		outfile.write('\t' + str(linecount))
		for key in keys:
			outfile.write('\t' + str(float(groupdict[group][key])/annotateddict[group][key]))
		outfile.write('\n')

