import argparse
from argparse import RawDescriptionHelpFormatter
import os
import os.path
import time
import subprocess
import signal

def getOptions():
        """Get arguments"""
        description="""This script can be used to get arg versus organism capturing orgs on the 
        contig and contigs up to N links away via Hi-C reads"""
        parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-i','--input',dest='inhandle',action='store',required=True,type=str, help='Kraken report', metavar='KRAKEN')
        parser.add_argument('-o','--output',dest='outhandle',action='store',required=True,type=str, help='Best hit with taxid for the bin', metavar='BESTHIT')
        args = parser.parse_args()
        return(args)

args = getOptions()
inhandle = args.inhandle
outhandle = args.outhandle
rankdict = {"species": 7, "genus":6, "family":5, "order":4,"class":3,"phylum":2,"superkingdom":1,"no rank":0, "unclassified": -1}

markerlist = ['k','p','c','o','f','g','s']
#|-__root|-__cellular_organisms|k__Bacteria|-__Terrabacteria_group|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Anaerotignum|s__Anaerotignum_propionicum
def parse_taxonomy(taxonomy):
	base_marker = {'k':'','p':'','c':'','o':'','f':'','g':'','s':''}
	taxonomy_fields = taxonomy.split('|')
	for field in taxonomy_fields[1:]:
		if field != '':
			marker = field.split('__')[0]
			taxa = field.split('__')[1]
			if marker in markerlist:
				base_marker[marker] = taxa
	parsedtaxonomy = ''
	for marker in markerlist:
		parsedtaxonomy = parsedtaxonomy + marker + '__' + base_marker[marker]
		if marker != 's':
			parsedtaxonomy = parsedtaxonomy + '; '
		else:
			parsedtaxonomy = parsedtaxonomy + ';'
	return parsedtaxonomy
			



#besttuple (percent,taxid, rank, ranknum) 
best = (0, '0', 'unclassified', -1)
unclassified = 0
with open(inhandle) as infile, open(outhandle+'.contamination','w') as outfile:
	for line in infile:
		if ((line[0] != '#') and (line[0] !='%') and (line[0] !='\n')):
			pct = float(line.split('\t')[0])
			taxid = line.split('\t')[6]
			rank = line.split('\t')[7]
			taxa = line.strip().split('\t')[8]
			#figure out unclassified
			if taxa == 'unclassified':
				unclassified = pct
			if pct>50:
				if rank in rankdict.keys():
					ranknum = rankdict[rank]
					if ranknum > best[3]:
						best = (pct, taxid, rank, ranknum)
	besttaxid = best[1]
	taxonomy=outhandle + '.krakentaxa'
	#rerunning so these should already be computed! hopefully
	cmd = '/workdir/users/agk85/tools/krakenuniq/query_taxdb /workdir/users/agk85/tools/krakenuniq/DB_subset/taxDB ' +  besttaxid  + ' > ' + taxonomy
	p = subprocess.Popen(cmd, shell=True,preexec_fn=os.setpgrp)
	while not (os.path.exists(taxonomy) and os.path.getsize(taxonomy) > 0):
		time.sleep(1)	
	with open(taxonomy, 'r') as taxfile:
		line = taxfile.readline()
		kraken_taxonomy = line.strip().split('\t')[1]
		parsed_taxonomy = parse_taxonomy(kraken_taxonomy)
	print(inhandle)
	print(parsed_taxonomy)
	contamination = str(abs(round(100.0-best[0]-unclassified,2)))
	binname = inhandle.split('.fa')[0].split('/')[-1]
	outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(binname, str(best[0]), best[1], best[2],best[3],contamination, parsed_taxonomy))

#os.kill(os.getpgid(p.pid), signal.SIGTERM)
