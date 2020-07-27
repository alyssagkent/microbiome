#USAGE python script /workdir/users/agk85/CDC
#goal is to create cleanliness plots r1 vs r2 given the taxonomic information for each contig they map to
#does it with mges and without mges

import numpy as np
import matplotlib.pyplot as plt
import glob
import collections
import sys
from Bio import SeqIO

#tasks
#1. 
#2.
#3.
#4.
base = sys.argv[1]
pidthresh = '0'

def best_org(inline):
	besttaxa = inline.split('\t')[-1].strip()
	return besttaxa



tbldict = {}
tblpaths = glob.glob(base + '/combo_tables/metagenomes/U*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()




hicdict = collections.defaultdict(dict)
hicdict_mge = collections.defaultdict(dict)
organisms = []
organisms2 = []
hicpaths = glob.glob(base + '/hicpro/output/*/hic_results/data/*/U*_trans_primary_' + pidthresh + '_noeuks.txt')
for pathway in hicpaths:
	with open(pathway) as hic:
		for line in hic:
			scf1 = line.split('\t')[2]
			scf2 = line.split('\t')[8]
			annot1 = tbldict[scf1]
			annot2 = tbldict[scf2]
			if((best_org(annot1)!='.') and (best_org(annot2)!='.')):
				if(annot1.split('\t')[-2] !='mge' and annot2.split('\t')[-2] !='mge'):
					taxa1 = best_org(annot1)
					taxa2 = best_org(annot2)	
					organisms.append(taxa1)
					organisms.append(taxa2)
					try:
						hicdict[taxa1].append(taxa2)
					except:
						hicdict[taxa1] = [taxa2]
				else:
					taxa1 = best_org(annot1)
					taxa2 = best_org(annot2)	
					organisms2.append(taxa1)
					organisms2.append(taxa2)
					try:
						hicdict_mge[taxa1].append(taxa2)
					except:
						hicdict_mge[taxa1] = [taxa2]
			

orglist = list(set(organisms))
orglist.sort()
header = ''
for org in orglist:
	header = header +'\t' + org

header = header + '\n'


delims = ['k__','p__','c__','o__','f__','g__','s__']
delim_names=['kingdom','phylum','class','order','family','genus','species']

def conflicted(levels):
	taxa  = list(set(levels))
	if '' in taxa:
		taxa.remove('')
	if len(taxa)>1:
		conflict = 1
	else:
		conflict = 0
	return conflict


def levelconflict(taxalist):
	flag = 0
	"""goal is to tell if at any level there were multiple non-empty taxa converted to empty"""
	for delim in delims:
		levs = list(set([taxonomy.split(delim)[1].split(';')[0] for taxonomy in taxalist if delim in taxonomy]))
		if conflicted(levs):
			flag = 1
	return flag 

mismatched = 0
total = 0

mismatched_levels = {}
total_levels= {}
for delim in delims:
	mismatched_levels[delim] = 0
	total_levels[delim] = 0

hicpaths = glob.glob(base + '/hicpro/output/*/hic_results/data/*/U*_trans_primary_' + pidthresh + '_noeuks.txt')
for pathway in hicpaths:
        with open(pathway) as hic:
                for line in hic:
                        scf1 = line.split('\t')[2]
                        scf2 = line.split('\t')[8]
                        annot1 = tbldict[scf1]
                        annot2 = tbldict[scf2]
                        if((best_org(annot1)!='.') and (best_org(annot2)!='.')):
				if(1==1):
                                #if(annot1.split('\t')[-2] !='mge' and annot2.split('\t')[-2] !='mge'):
                                        taxa1 = best_org(annot1)
                                        taxa2 = best_org(annot2)
                                        if  levelconflict([taxa1, taxa2]):
                                                mismatched += 1
						total += 1
					else:
						total += 1
					for delim in delims:
						t1 = taxa1.split(delim)[1].split(';')[0]
						t2 = taxa2.split(delim)[1].split(';')[0]
						if (t1!='' and t2!= ''):
							total_levels[delim] += 1
							if conflicted([t1,t2]):
								mismatched_levels[delim] += 1
						
print("mismatched")
print(mismatched)
print("out of")
print(total)


for delim in delims:
	print("mismatch",delim, mismatched_levels[delim])
	print("total",delim, total_levels[delim])




mismatch = {}
totcount = {}
for delim in delims:
	mismatch[delim] = 0
	totcount[delim] = 0


count = 0
#you might have a rule break in some case where you don't get every org in the taxa 1 position
outhandle = base + '/hicpro/figures/U_org_v_org_' + pidthresh + '.txt'
mismatchhandle = base + '/hicpro/figures/U_org_v_org_' + pidthresh + '_mismatches_variablecount.txt'
with open(outhandle,'w') as outfile:
	with open(mismatchhandle, 'w') as mismatchfile:
		outfile.write(header)
		for org in orglist:
			orgmismatch = 0
			orgcount = 0
			outfile.write(org)
			try:
				h = hicdict[org]
				for org2 in orglist:
					outfile.write('\t' + str(h.count(org2)))
					orgcount = orgcount + h.count(org2)
					count = count + h.count(org2)
					if org == org2:
						for delim in delims:
							totcount[delim] = totcount[delim] + h.count(org2)
					if org != org2:
						for delim in delims:
							t1 = org.split(delim)[0] + delim + org.split(delim)[1].split(';')[0] + ';'
							t2 = org2.split(delim)[0] + delim + org2.split(delim)[1].split(';')[0] + ';'
							if levelconflict([t1, t2]):
								mismatch[delim] = mismatch[delim] + h.count(org2)
								totcount[delim] = totcount[delim] + h.count(org2)
								if delim =='p__' and h.count(org2)>0:
									print org
									print org2
									print(str(h.count(org2)))
						if levelconflict([org, org2]):
							orgmismatch = orgmismatch + h.count(org2)
			except AttributeError:
				for org2 in orglist:
					outfile.write('\t0')
			outfile.write('\n')
			try:
				div = orgmismatch/float(orgcount)
			except ZeroDivisionError:
				div = 0
			mismatchfile.write(org + '\t' + str(orgmismatch) + '\t' + str(orgcount) + '\t' + str(div) + '\n')
		for delim in delims:
			try:
				div = mismatch[delim]/float(count)
			except ZeroDivisionError:
				div = 0
			mismatchfile.write(delim + '\t' + str(mismatch[delim]) + '\t' + str(totcount[delim]) + '\t' + str(div) + '\n')


orglist = list(set(organisms2))
orglist.sort()
header = ''
for org in orglist:
	header = header +'\t' + org

header = header + '\n'
outhandle = base + '/hicpro/figures/U_org_v_org_' + pidthresh + '_mge.txt'
with open(outhandle,'w') as outfile:
	outfile.write(header)
	for org1 in orglist:
		outfile.write(org1)
		try:
			h = hicdict_mge[org1]
			for org2 in orglist:
				outfile.write('\t' + str(h.count(org2)))
		except AttributeError:
			for org2 in orglist:
				outfile.write('\t0')
		outfile.write('\n')
	
	

