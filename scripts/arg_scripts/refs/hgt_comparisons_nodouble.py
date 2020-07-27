#how to make file
import sys
inhandle = sys.argv[1]
outhandle = sys.argv[2]

hgtdict = {}
patients = ['B314','B316','B320','B331','B335','B357','B370']
levels = ['genus','family','order','class','phylum','kingdom']
for level in levels:
	hgtdict[level] = {}
	for pat in patients:
		hgtdict[level][pat] = {}


with open(inhandle,'r') as infile:
	for line in infile:
		fields = line.strip().split('\t')
		pat1 = fields[0]
		pat2 = fields[1]
		taxalevel = fields[5]
		hgtdict[taxalevel][pat1][pat2] = line


finished = []
with open(outhandle,'w') as outfile:
		for level in hgtdict:
			leveldict = hgtdict[level]
			for pat1 in leveldict:
				pat1dict = leveldict[pat1] 
				for pat2 in pat1dict:
					if level + pat2 + pat1 not in finished:
						outfile.write(hgtdict[level][pat1][pat2])
						finished.append(level+pat1+pat2)
