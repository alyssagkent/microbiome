

#percent mobile gene patient overlap?
#goal run through mobile genes, and figure out what % have more than one patient, more than 2 more than 3?


desc = {}
with open('mge_95_nr.fna.clstr.tbl.annot') as infile:
	for line in infile:
		cluster =  line.split('\t')[0]
		desc[cluster] = line.split('\t')[1]





with open ('mge_95_nr.fna.clstr.tbl') as infile:
	mt0 = 0
	mt1 = 0
	mt2 = 0
	mt3 = 0
	mt4 = 0
	mt5 = 0
	mt6 = 0
	for line in infile:
		cluster = line.split('\t')[0]
		genes = line.split('\t')[2].split(',')
		pats = []
		for gene in genes:
			pats.append(gene.split('-')[0])
		patset = len(list(set(pats)))
		mt0 += 1
		if patset >1:
			mt1 += 1
			if patset >2:
				mt2 += 1
				if patset>3:
					mt3 += 1
					if patset>4:
						mt4 += 1
						if patset>5:
							mt5 += 1
							if patset>6:
								mt6 += 1
								print(desc[cluster])


with open('mobile_gene_patient_overlap.txt','w') as outfile:
	outfile.write('1+ gene clusters:{0}\t{1}\n'.format(str(mt0), str(round(100*float(mt0)/mt0,4))))
	outfile.write('2+ gene clusters:{0}\t{1}\n'.format(str(mt1), str(round(100*float(mt1)/mt0,4))))
	outfile.write('3+ gene clusters:{0}\t{1}\n'.format(str(mt2), str(round(100*float(mt2)/mt0,4))))
	outfile.write('4+ gene clusters:{0}\t{1}\n'.format(str(mt3), str(round(100*float(mt3)/mt0,4))))
	outfile.write('5+ gene clusters:{0}\t{1}\n'.format(str(mt4), str(round(100*float(mt4)/mt0,4))))
	outfile.write('6+ gene clusters:{0}\t{1}\n'.format(str(mt5), str(round(100*float(mt5)/mt0,4))))
	outfile.write('7+ gene clusters:{0}\t{1}\n'.format(str(mt6), str(round(100*float(mt6)/mt0,4))))
	
	
			
