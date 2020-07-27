import glob
from Bio import SeqIO
import collections
import argparse
from argparse import RawDescriptionHelpFormatter

def getOptions():
        """Get arguments"""
        description="""This script takes the taxa from the gaemr contigs, assigns to the scaffolds, reports the taxa with the longest assignment"""
        parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-s','--scaffold',dest='scaffold',action='store',required=True, type=str,help="scaffold file",metavar='SCAFFOLD')
	parser.add_argument('-n', '--name', dest='name', action='store', required=True, type=str, help='Sample Name [Required]', metavar='NAME') 
        parser.add_argument('-i', '--infolder', dest='infolder',action = 'store', required=True, type=str, help='Folder where gaemr top folders are /workdir/users/agk85/CDC/gaemr/metagenomes5/', metavar='INFOLDER')
	args = parser.parse_args()
        return(args)

args = getOptions()
name = args.name
infolder = args.infolder
ids = []
scffile = args.scaffold
with open(scffile,'r') as infile:
        for line in infile:
                if line[0] == '>':
                        id = line.strip().split(' ')[0].split('>')[1]
                        ids.append(id)
                        

def get_taxa(taxa, splitter):
	"""Define how to split the taxa"""
	try:
		taxon = taxa.split(splitter+'=')[1].split(';')[0]
	except IndexError:
		taxon = ''
	return taxon
	             

#for each contig get the taxonomy
contig_taxa = {}
contig_length = {}
taxonomies = infolder + name + '/gaemr/table/' + name + '.blast_hit_taxonomy.table.txt'
with open(taxonomies) as taxafile:
	header = taxafile.readline()
	header = taxafile.readline()
	for line in taxafile:
		fields = line.strip().split(' | ')
		contig = fields[0]
		length = fields[1]
		QueryHitLen = fields[2]
		contig_length[contig] = int(QueryHitLen.replace(',',''))
		pct_covered = fields[3]
		taxa = fields[4]
		if taxa != 'domain=NO_HIT':
			#domain=cellular organisms;superkingdom=Bacteria;domain=Terrabacteria group;phylum=Firmicutes;class=Negativicutes;order=Veillonellales;family=Veillonellaceae;genus=Dialister;species=Dialister sp. Marseille-P5638
			k = get_taxa(taxa, 'superkingdom')
			p = get_taxa(taxa, 'phylum')
			c = get_taxa(taxa, 'class')
			o = get_taxa(taxa, 'order')
			f = get_taxa(taxa, 'family')
			g = get_taxa(taxa, 'genus')
			s = get_taxa(taxa, 'species')
			if s != 'uncultured bacterium':
				taxon = 'k__' + k + '; p__' + p + '; c__' + c + '; o__' + o + '; f__' + f + '; g__' + g + '; s__' + s + ';'
				contig_taxa[contig] = taxon
			else:
				contig_taxa[contig] = '.'
		else:
			contig_taxa[contig] = '.'


#connect the contig id to the scaffold id
contig_scaffold = {}
contigscaffold = infolder + name + '/gaemr/table/' + name + '.contig_detail.table.txt'
with open(contigscaffold) as linkerfile:
	header = linkerfile.readline()
	header = linkerfile.readline()
	for line in linkerfile:
		contig = line.split(' | ')[0]
		scaffold = line.split(' | ')[1]
		pct = line.split(' | ')[6].strip()
		contig_scaffold[contig] = (scaffold, pct)

#connect the scaffold ids to the gaemr scaffold ids
scfs = infolder + name + '/gaemr/work/' + name + '.scaffolds.fasta'
scfdict = SeqIO.to_dict(SeqIO.parse(scfs, "fasta"))
scaffolddict = SeqIO.to_dict(SeqIO.parse(scffile, "fasta"))
scflink = {}
count = 0
for rec in SeqIO.parse(scfs, "fasta"):
	s1 = str(scfdict[rec.id].seq)
	s2 = str(scaffolddict[ids[count]].seq)
	scflink[rec.id] = ids[count]
	#if s1 != s2:
	#	print rec.id
	#	print ids[count]
	count = count + 1


#so keep track of length of each taxa with a length dictionary
scftaxalen = {}

for contig in contig_taxa.keys():
	#get the scaffold that contig is associated with
	gaemrscf = contig_scaffold[contig][0]
	gaemrpct = float(contig_scaffold[contig][1])
	#get the taxa that contig is associated with
	taxa = contig_taxa[contig]
	#get the length matched
	taxalength = contig_length[contig]
	taxahitlength = contig_length[contig]*gaemrpct/100
	try:
		try:
			scftaxalen[gaemrscf][taxa][0] = scftaxalen[gaemrscf][taxa][0] + taxahitlength
			scftaxalen[gaemrscf][taxa][1] = scftaxalen[gaemrscf][taxa][1] + taxalength
		except KeyError:
			scftaxalen[gaemrscf][taxa] = (taxahitlength, taxalength)
	except:
		scftaxalen[gaemrscf] = {}
		scftaxalen[gaemrscf][taxa] = (taxahitlength, taxalength)


#ok so now that you have that dictionary make a new one with a single taxa with the most
scftaxa = {}
for gaemrscf in scftaxalen.keys():
	longesttaxa = '.'
	l = 0
	for taxa in scftaxalen[gaemrscf].keys():
		newl = scftaxalen[gaemrscf][taxa][0]
		if newl >= l:
			longesttaxa = taxa
			l = newl
			total = scftaxalen[gaemrscf][taxa][1]
			if l != 0:
				percent = float(newl)/total	
			else:
				percent = 0
	scftaxa[gaemrscf] = (longesttaxa,percent)



outhandle = infolder + name + '/gaemr/table/' + name + '.scf.taxa.percent.txt'
with open(outhandle,'w') as outfile:
	for gaemrscf in scftaxa.keys():
		realscf = scflink[gaemrscf]
		taxa = scftaxa[gaemrscf][0]
		percent = scftaxa[gaemrscf][1]
		replace1 = taxa.replace(',','')
		replace2 = replace1.replace("'", "")
		replace3 = replace2.replace("#","")
		outfile.write(realscf + '\t' + replace3 + '\t' + str(percent) + '\n')

