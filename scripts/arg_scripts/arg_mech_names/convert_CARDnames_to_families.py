#I added three terms to this file (the first three) that were no longer in this version using their closest counterparts
import sys
import os

folder = sys.argv[1]
argpid =sys.argv[2]
folder='CDC2'
argpid='99'

#do this before
#os.system("cat ~/agk/" + folder + "/card/metagenomes/*/*.txt | grep -v '^ORF_ID' > ~/agk/" + folder + "/arg_v_org/metagenomes/all_card.txt")

names = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/mapping/bwa_alignments_' + argpid + '_' + argpid + '/arg_v_samp_' + argpid + '_' + argpid + '_names.txt'
mechhandle = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/arg_v_samp_' + argpid + '_' + argpid + '_names_mech.txt'
clstrhandle = '/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/args_' + argpid + '_nr.fna.clstr'


#you have these already 
#resfam_handle = '/workdir/users/agk85/press/scripts/arg_scripts/arg_mech_names/all_hic_resfams.txt'
#card_handle = '/workdir/users/agk85/press/scripts/arg_scripts/arg_mech_names/all_hic_card.txt'

resfam_handle='/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/resfams_prot_id_type.txt'
card_handle='/workdir/users/agk85/' + folder + '/arg_v_org/metagenomes/all_card.txt'

#this is the intermediate outhandle---just leave them here
resfam_mechhandle = '/workdir/users/agk85/CDC2/scripts/arg_scripts/arg_mech_names/all_hic_resfams_mech.txt'
card_mechhandle = '/workdir/users/agk85/CDC2/scripts/arg_scripts/arg_mech_names/all_hic_card_mech.txt'


obofile ='/workdir/users/agk85/CDC2/scripts/arg_scripts/arg_mech_names/aro.obo.agk'
resfams_metadata = '/workdir/users/agk85/CDC2/scripts/arg_scripts/arg_mech_names/180102_resfams_metadata_updated_v1.2.2.txt'
megares_categories = '/workdir/users/agk85/CDC2/scripts/arg_scripts/arg_mech_names/megares_annotations_v1.01.csv'
megares_maphandle = '/workdir/users/agk85/CDC2/scripts/arg_scripts/arg_mech_names/megares_to_external_header_mappings_v1.01.tsv'



def get_stanzas(obofile):
	stanzas = []
	with open(obofile) as obo:
		stanza = ''
		for line in obo:
			if line != '\n':
					stanza = stanza + line
			if line == '\n':
				stanzas.append(stanza)
				stanza = ''
	return stanzas

arodict = {}
aroisadict = {}
stanzas = get_stanzas(obofile)
for stanza in stanzas:
	if '[Term]' in stanza:
		aroid = stanza.split('\nid: ')[1].split('\n')[0]
		arodict[aroid] = []
		aroisadict[aroid] = []
		if 'is_a' in stanza:
			aroisas = stanza.split('\nis_a: ')
			for links in aroisas[1:]:
				aroisa = links.split(' ! ')[0]
				aromech = links.split(' ! ')[1].split('\n')[0]
				arodict[aroid].append(aroisa)
				aroisadict[aroid].append(aromech)

#replace efflux
replacements = {"RND Antibiotic Efflux":"Efflux",
"Multi-drug efflux pumps":"Efflux",
"ABC Transporter":"Efflux",
"Other Efflux":"Efflux",
"MFS Transporter":"Efflux",
"Gene Modulating Resistance,Gylcopeptide Resistance":"Gylcopeptide Resistance",
"Gene Modulating Resistance,Efflux":"Efflux",
"Lipid A modification":"Phosphoethanolamine Transferase",
"Fluoroquinolone-resistant DNA topoisomerases":"Quinolone Resistance",
"Macrolide phosphotransferase":"Phosphotransferase",
"Macrolide Phosphotransferase":"Phosphotransferase",
"Macrolide phosphotransferases":"Phosphotransferase",
"MDR regulator":"Gene Modulating Resistance",
"Phenicol 23S rRNA methyltransferases":"rRNA Methyltransferases",
"Aminoglycoside N-acetyltransferases":"Acetyltransferase",
"Daptomycin-resistant liaFSR":"Daptomycin-Resistant liaFSR",
"Dihydrofolate reductase":"Dihydrofolate Reductase",
"Phosphoethanolamine transferase":"Phosphoethanolamine transferase",
"Sulfonamide-resistant dihydropteroate synthases":"Sulfonamide Resistant Dihydropteroate Synthases",
"Undecaprenyl pyrophosphate phosphatase":"Undecaprenyl Pyrophosphate Phosphatase",
"Fosfomycin thiol transferases":"Fosfomycin Thiol Transferase",
"Fosfomycin transporter":"Fosfomycin Transporter",
"Fosfomycin Thiol Transferases":"Fosfomycin Thiol Transferase",
"Isoleucyl-tRNA transferase":"Isoleucyl-tRNA Transferase",
"Lincosamide nucleotidyltransferases":"Lincosamide Nucleotidyltransferases",
"Efflux,Tetracycline MFS Efflux":"Tetracycline MFS Efflux",
"Lincosamide Nucleotidyltransferases":"Nucleotidyltransferase",
"rRNA Methyltransferases": "rRNA Methyltransferase",
"Fosfomycin Thiol Transferases":"Fosfomycin Thiol Transferase",
"VanG-type resistance protein":"Gylcopeptide Resistance",
"Beta-Lactamase-Class A,Gene Modulating Resistance":"Beta-Lactamase-Class A"
}

#resfams mechanisms
resfam_dict = {}
with open(resfams_metadata) as r:
	header = r.readline()
	for line in r:
		ResfamID=line.split('\t')[0]
		name = line.split('\t')[1]
		mech = line.split('\t')[6]
		tetextrainfo = line.split('\t')[8].strip()
		if mech == 'Beta-Lactamase':
			mech = mech + '-' + line.split('\t')[7]
		if tetextrainfo != '':
			mech = line.split('\t')[8].strip()
		if mech in replacements.keys():
			mech = replacements[mech]
		resfam_dict[name] = mech

resfam_geneids_mech = {}
with open(resfam_mechhandle,'w') as resoutfile:
	with open(resfam_handle) as resfam_map:
		for line in resfam_map:
			geneid = line.split('\t')[0]
			resname = line.strip().split('\t')[1]
			mech = resfam_dict[resname]
			resoutfile.write(geneid + '\t' + resname + '\t' + mech + '\n')
			if mech in replacements.keys():
				mech = replacements[mech]
			resfam_geneids_mech[geneid]= mech


resfam_card_dict = {}
with open(resfams_metadata) as r:
	header = r.readline()
	for line in r:
		aros = line.split('\t')[3].split('; ')
		mech = line.split('\t')[6]
		tetextrainfo = line.split('\t')[8].strip()
		if mech == 'Beta-Lactamase':
			mech = mech + '-' + line.split('\t')[7]
		if tetextrainfo != '':
			mech = line.split('\t')[8].strip()
		for aro in aros:
			if mech in replacements.keys():
				mech = replacements[mech]
			resfam_card_dict[aro] = mech


#megares category conversions
megares_mech = {}
with open(megares_categories) as infile:
	header = infile.readline()
	for line in infile:
		megares = line.split(',')[0]
		mech = line.split(',')[2]
		if mech in replacements.keys():
			mech = replacements[mech]
		megares_mech[megares] = mech


arotomech = {}
arotoclass = {}
with open(megares_maphandle) as reffile:
	header = reffile.readline()
	for line in reffile:
		source = line.split('\t')[0]
		if source == 'CARD':
			megares = line.split('\t')[1]
			cardaroid = line.strip().split('\t')[2].split('|')[3]
			arotomech[cardaroid] = megares_mech[megares]


card_exceptions = {'ARO:3004281':'General Bacterial Porin Conferring Beta-lactam Resistance','ARO:3000446':'Isoleucyl-tRNA Transferase','ARO:3004248':'Fosfomycin Transporter','ARO:3004268':'Phosphoethanolamine Transferase','ARO:3000333':'Macrolide Phosphotransferase','ARO:3004247':'Fosfomycin Transporter','ARO:3002866':"Dihydrofolate Reductase",'ARO:3003210':"Fosfomycin Thiol Transferase",'ARO:3004045':"Determinant of Isoniazid Resistance"}

	
count = 0
dealwith = []
with open(card_handle) as cardfile:
	with open(card_mechhandle, 'w') as cardoutfile:	
		cardoutfile.write('gene\tname\tbest_mech\tresfam_mech\tmegares_mech\tcard_mech\tcardisamech\n')
		for line in cardfile:
			annotations = {'resfam':[],'megares':[],'card':[], 'cardisa':[]}
			gene = line.split('\t')[0].split(' ')[0]
			if gene in resfam_geneids_mech.keys():
				annotations['resfam'].append(resfam_geneids_mech[gene])
			print(line.split('\t')[8])
			name = line.split('\t')[8]
			aros = line.split('\t')[10].split(', ')
			alllinkedaros = []
			for aro in aros:
				alllinkedaros.append(aro)
				try:
					morearos = arodict[aro]
					newmorearos = []
					for newaro in morearos:
						newmorearos = arodict[newaro]
					alllinkedaros = alllinkedaros + morearos + newmorearos
				except KeyError:
					a = 1
			for aro in alllinkedaros:
				if aro in resfam_card_dict.keys():
					annotations['resfam'].append(resfam_card_dict[aro])
				if aro in arotomech.keys():
						annotations['megares'].append(arotomech[aro])
				if aro in card_exceptions.keys():
					annotations['card'].append(card_exceptions[aro])
				if aro in aroisadict.keys():
					for item in aroisadict[aro]:
						annotations['cardisa'].append(item)
			resfammech = ','.join(list(set(annotations['resfam'])))
			megaresmech = ','.join(list(set(annotations['megares'])))
			cardmech = ','.join(list(set(annotations['card'])))
			cardisamech = ','.join(list(set(annotations['cardisa'])))
			#get the best mech resfam then megares then card in terms of importance
			bestmech = ''
			if cardmech != '':
				bestmech = cardmech
			if megaresmech != '':
				bestmech = megaresmech
			if resfammech !='':
				bestmech = resfammech
			if resfammech+megaresmech+cardmech == '':
				count = count + 1
				print(name)
				print(alllinkedaros)
				print('')
				dealwith.append((name,aro))
			if bestmech in replacements.keys():
				bestmech = replacements[bestmech]
			cardoutfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(gene, name, bestmech, resfammech,megaresmech, cardmech, cardisamech))

todealwith = set(dealwith)
todealwith


####################cluster to 
numcluster_map = {}
with open(clstrhandle) as infile:
	for line in infile:
		if line[0] == '>':
			cluster = line.strip().split(' ')[1]
		else:
			prot =  line.split('>')[1].split('...')[0]
			if '*' in line:	
				numcluster_map[prot] = cluster


resargmech = {}
with open(resfam_mechhandle) as resfile:
	for line in resfile:
		gene = line.split('\t')[0]
		argmech = line.split('\t')[2].strip()
		resargmech[gene] = argmech


cardargmech = {}	
with open(card_mechhandle) as cardfile:
	header = cardfile.readline()
	for line in cardfile:
		gene = line.split('\t')[0].split(' ')[0]
		argmech = line.split('\t')[2]
		cardargmech[gene] = argmech


subcluster = {"Beta-Lactamase-Class A":"Beta-Lactamase",
"Beta-Lactamase-Class B":"Beta-Lactamase",
"Beta-Lactamase-Class C":"Beta-Lactamase",
"Beta-Lactamase-Class D":"Beta-Lactamase",
"Penicillin binding protein":"Beta-Lactam Resistance",
"General Bacterial Porin Conferring Beta-lactam Resistance":"Beta-Lactam Resistance",
"Tetracycline MFS Efflux":"Efflux",
"Tetracycline Inactivation":"Tetracycline Resistance",
"Tetracycline Ribosomal Protection":"Tetracycline Resistance",
"Fosfomycin Thiol Transferases":"Fosfomycin Resistance",
"Fosfomycin Transporter":"Fosfomycin Resistance",
"Efflux":"Efflux",
"Acetyltransferase":"Acetyltransferase",
"Nucleotidyltransferase":"Nucleotidyltransferase",
"Phosphotransferase":"Phosphotransferase",
"Phosphoethanolamine Transferase":"Phosphoethanolamine Transferase",
"rRNA Methyltransferase":"rRNA Methyltransferase",
"Gene Modulating Resistance":"Gene Modulating Resistance",
"Gylcopeptide Resistance":"Gylcopeptide Resistance",
"Quinolone Resistance":"Quinolone Resistance"
}

with open(names) as infile:
	with open(mechhandle,'w') as outfile:
		header = infile.readline()
		outfile.write('Cluster' + '\t' + header.strip() + '\t' + 'Mechanism\tSub_mechanism\n')
		for line in infile:
			gene = line.split('\t')[0]
			cluster = numcluster_map[gene]
			if gene in resargmech.keys():
				bestmech = resargmech[gene]
			else:
				bestmech = cardargmech[gene]
			if bestmech in subcluster.keys():
				submech = subcluster[bestmech]
			else:
				submech = "Other"
			outfile.write(cluster + '\t' + line.strip() + '\t' + bestmech + '\t' + submech+'\n')

