#goal is to produce a arg versus organism table with extra info about the number of genes, samples, patients that 
#USAGE python script /workdir/users/agk85/CDC 90
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
pidthresh = sys.argv[2]
#pidthresh = '90'

def best_org(inline):
	besttaxa = inline.split('\t')[-1].strip()
	return besttaxa



tbldict = {}
tblpaths = glob.glob(base + '/combo_tables/metagenomes4/*_master_scf_table.txt')
for tblfile in tblpaths:
	with open(tblfile) as t:
		header = t.readline()
		for line in t:
			tbldict[line.split('\t')[0]] = line.strip()

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

hicpaths = glob.glob(base + '/newhic/mapping/*/*_trans_primary_' + pidthresh + '_noeuks.txt')
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
                                #if(annot1.split('\t')[-2] =='mge' and annot2.split('\t')[-2] =='mge'):
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


print("HIC_level_mismatch_total","any",mismatched, total)
for delim in delims:
        print("HIC_level_mismatch_total",delim, mismatched_levels[delim], total_levels[delim])
						
mismatched = 0
total = 0
mismatched_levels = {}
total_levels= {}
for delim in delims:
	mismatched_levels[delim] = 0
	total_levels[delim] = 0

hicpaths = glob.glob(base + '/newhic/mapping_mgm/*/*_trans_primary_' + pidthresh + '_noeuks.txt')
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
				#if(annot1.split('\t')[-2] =='mge' and annot2.split('\t')[-2] =='mge'):
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

print("MGM_level_mismatch_total","any",mismatched, total)
for delim in delims:
        print("MGM_level_mismatch_total",delim, mismatched_levels[delim], total_levels[delim])
