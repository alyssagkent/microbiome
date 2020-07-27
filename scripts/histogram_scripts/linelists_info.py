inhandle = 'Together_das_2_mgetaxa_linelist.txt'
outhandle = inhandle.split('linelist')[0] + 'linelist_info.txt'

mgedict = {}
with open(reffile) as reffile:
	for line in reffile:
		mgeid = line.split('\t')[0]
		mgetype = line.split('\t')[3] #this is any gene in the group (lots of multis)
		mgedict[mgeid] = mgetype


with open(inhandle) as infile, open(outhandle,'w') as outfile:
	for line in infile:
		mgeid = line.split('\t')[6]
		mgetype = mgedict[mgeid]
		outfile.write(line.strip() + '\t' + mgetype + '\n')
