#quick script to get me all of the connections in a patient
import glob

reffile = '/workdir/users/agk85/CDC2/mobile/metagenomes/mge_99_nr.fna.clstr.tbl.annot'
mgedict = {}
with open(reffile) as reffile:
        for line in reffile:
                mgeid = line.split('\t')[0]
                mgetype = line.split('\t')[3] #this is any gene in the group (lots of multis)
                mgedict[mgeid] = mgetype

reffile = '/workdir/users/agk85/CDC2/args/arg_v_samp_99_99_names_mech.txt'
argdict = {}
with open(reffile) as reffile:
	header = reffile.readline()
	for line in reffile:
		argid = line.split('\t')[0]
		argtype = line.strip().split('\t')[6]
		argdict[argid] = argtype

infiles = glob.glob('/workdir/users/agk85/CDC2/bins/*_das_*_*taxa.txt')
for inhandle in infiles:
	if 'all' not in inhandle:
		print(inhandle)
		if 'arg' in inhandle:
			genedict = argdict
		if 'mge' in inhandle:
			genedict = mgedict
		outhandle = inhandle.split('taxa.txt')[0] + 'taxa_linelist.txt'
		with open(inhandle) as infile, open(outhandle,'w') as outfile:
			header = infile.readline()
			#outfile.write('Min_contacts\tPatient\tThresh\tBintype\tResidency\tLevel\tGene\tGenetype\tTaxon\n')
			for line in infile:
				Count,Patient,Thresh,Bintype,Residency,Level,Gene,TaxaCount,Total_mges,Taxaset,Min_contacts = line.strip().split('\t')
				if (Level != 'k__') and (Thresh == 'contacts') and (Bintype == 'anybin'):
					taxa = Taxaset.split(',')
					genetype = genedict[Gene]
					uniquetaxa = list(set(taxa))
					for taxon in uniquetaxa:
						outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(Min_contacts,Patient,Thresh,Bintype,Residency,Level,Gene,genetype,taxon))
			
	


#rm all*linelist.txt
#cat *_das_2_argtaxa_linelist.txt > Together_das_2_argtaxa_linelist.txt
#cat *_das_5_argtaxa_linelist.txt > Together_das_5_argtaxa_linelist.txt
#cat *_das_2_mgetaxa_linelist.txt > Together_das_2_mgetaxa_linelist.txt
#cat *_das_5_mgetaxa_linelist.txt > Together_das_5_mgetaxa_linelist.txt

