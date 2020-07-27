card_dict = {}
with open('/workdir/users/agk85/CDC2/args/all_card.txt') as infile:
	for line in infile:
		name = line.split('\t')[8]
		#names = line.split('\t')[11].split(', ')
		aros = line.split('\t')[10].split(', ')
		card_dict[name] = aros
			#outfile.write('{0}\t{1}\t{2}\n'.format(gene,name,aro))

resfams_dict = {}
with open('/workdir/users/agk85/CDC2/args/resfams_mapping_file.txt') as reffile:
	for line in reffile:
		name = line.split('\t')[1]
		aros = line.strip().split('\t')[2].split('; ')
		resfams_dict[name] = aros

arodrug = {}
with open('/workdir/users/agk85/CDC2/args/aro_name_drug.txt') as dfile:
	for line in dfile:
		aro = line.split('\t')[0]
		name = line.split('\t')[1]
		drug = line.split('\t')[2]
		try:
			arodrug[aro].append(drug)
		except KeyError:
			arodrug[aro]= [drug]


with open('/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech_resfinder.txt') as infile, open('/workdir/users/agk85/CDC2/args/abundances/arg_v_samp_99_99_names_mech_refinder_drug.txt','w') as outfile:
	header = 'Cluster\tProtein\tName\tCARD\tResfams\tMechanism\tSub_mechanism\tResfinder_name\tResfinder_mech\tAROID\tDrug\n'
	outfile.write(header)
	for line in infile:
		cluster,protein,name,card,resfams,mechanism,sub_mechanism,resfinder_name,refinder_mech = line.strip().split('\t')
		flag = 0
		#try card
		if card != 'NA':
			try:
				card_aros = card_dict[card]
				for aro in card_aros:
					try:
						drugs= arodrug[aro]
						for drug in drugs:
							outfile.write(line.strip() + '\t'+aro + '\t' + drug + '\n')
						flag = 1
					except KeyError:
						a = 1
			except KeyError:
				a = 1
		#try resfams
		if resfams != 'NA':
			try:
				resfams_aros = resfams_dict[resfams]
				for aro in resfams_aros:
					try:
						drugs = arodrug[aro]
						for drug in drugs:
							outfile.write(line.strip() + '\t'+aro + '\t' + drug + '\n')
						flag = 1
					except KeyError:
						a = 1
			except KeyError:
				a=1
		#if neither
		if flag == 0:
			aro = 'NA'
			drug = 'NA'
			outfile.write(line.strip() + '\t'+aro + '\t' + drug + '\n')

