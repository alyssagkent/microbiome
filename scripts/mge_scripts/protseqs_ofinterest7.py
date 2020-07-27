#awk '$6 ~ /plasmid_/' B309-1_master_protein_table.txt
#
#USAGE: python protseqs_ofinterest.py
####################################################################################

from Bio import SeqIO
import os
import sys
import collections
name = sys.argv[1]
folder = sys.argv[2]

protdict = collections.defaultdict(dict)
mgedict = collections.defaultdict(dict)

def add_mge(mgedict, prot, mgetype, seq):
	mgedict[prot]['seq'] = seq
	try:
		mgedict[prot]['type'].append(mgetype)
	except KeyError:
		mgedict[prot]['type'] = [mgetype]


def update_with_scaffold(start, stop, scf, mgetype):
	with open('/workdir/users/agk85/'+folder+'/prodigal_excise/metagenomes/' + name + '/' + name + '_proteins.fna') as pf:
		for protline in pf:
			if protline[0] == '>':
				#check to make sure the scaffold is the same
				if scf.split('_')[2] == protline.split('_')[2]:
					#grab the beginning and end positions of the protein and the protein id
					b = int(protline.split(' # ')[1])
					e = int(protline.split(' # ')[2])
					prot = protline.split(' ')[0].split('>')[1]
					#if they overlap and if the overlapping section covers more than 50% of IS element length
					if (start < e and stop > b and (float(min(e-start, stop -b, stop-start, e-b)))/(e-b)>0.5):
						seq = protdict[prot]
						protid = protline.split('>')[1].split(' ')[0]
						add_mge(mgedict, protid, mgetype, seq)

print(name)
protfile = '/workdir/users/agk85/'+folder+'/prodigal_excise/metagenomes/' + name + '/' + name + '_proteins.fna'
for rec in SeqIO.parse(protfile, "fasta"):
	protdict[rec.id] = str(rec.seq)

#protdict_old = {}
#protfile = '/workdir/users/agk85/'+folder+'/prodigal/metagenomes/' + name + '/' + name + '_proteins.fna'
#for rec in SeqIO.parse(protfile, "fasta"):
#	protdict_old[rec.id] = str(rec.seq)	


#Relaxase
print('relaxase')
with open('/workdir/users/agk85/'+folder+'/plasmids/metagenomes/' + name + '/' + name + '_relaxase.txt','r') as relaxase:
	#relaxase_dict = {}
	for line in relaxase:
		prot = line.split('\t')[0]
		seq = protdict[prot]
		mgetype = 'plasmid'
		add_mge(mgedict, prot, mgetype, seq)
#phage--vogs
#print('vogs')
#with open('/workdir/users/agk85/'+folder+'/phage/metagenomes/' + name + '/' + name + '_vogs.txt.best.filter','r') as vogs:
#	for line in vogs:
#		prot = line.split('\t')[0]
#		seq = protdict[prot]
#		mgetype = 'phage'
#		add_mge(mgedict, prot, mgetype, seq)

#phage--phaster
print('phaster')
with open('/workdir/users/agk85/'+folder+'/phage/metagenomes/' + name + '/' + name + '_phaster_filter.out','r') as phaster:
	for line in phaster:
		prot = line.split('\t')[0]
		seq = protdict[prot]
		mgetype = 'phage'
		add_mge(mgedict, prot, mgetype, seq)

#aclame--all
#print('aclame')
##vir = phage, plasmid = plasmid, proph = phage
with open('/workdir/users/agk85/'+folder+'/aclame/metagenomes/' + name + '/' + name + '_aclame_filter.out','r') as aclame:
	for line in aclame:
		prot = line.split('\t')[0]
		elementtype = line.split('\t')[1].split(':')[1]
		seq = protdict[prot]
		if elementtype == 'plasmid':
			mgetype = 'plasmid'
			add_mge(mgedict, prot, mgetype, seq)
		if elementtype == 'proph':
			mgetype = 'phage'
			add_mge(mgedict, prot, mgetype, seq)
		if elementtype == 'vir':
			mgetype = 'phage'
			add_mge(mgedict, prot, mgetype, seq)

#phagefinder
print('phagefinder')
with open('/workdir/users/agk85/'+folder+'/prodigal_excise/metagenomes/' + name + '/' + name + '_phagefinder.txt','r') as phagefinder:
	for line in phagefinder:
		start = int(line.split('\t')[2])
		stop = int(line.strip().split('\t')[3])
		scf = line.split('\t')[0].split('>')[1]
		mgetype = 'phage'
		update_with_scaffold(start, stop, scf, mgetype)


print('isescan')
with open('/workdir/users/agk85/'+folder+'/iselements/metagenomes/prediction/' + name + '_scaffold.fasta.gff','r') as isescan:
	for line in isescan:
		if 'insertion_sequence' in line:
			start = int(line.split('\t')[3])
			stop = int(line.split('\t')[4])
			scf = line.split('\t')[0]
			mgetype = 'transposon'
			update_with_scaffold(start, stop, scf, mgetype)

print('plasmid_finder')
with open('/workdir/users/agk85/'+folder+'/plasmids/metagenomes/' + name + '/' + name + '_plasmidgenes_filter.out') as plasmid_finder: 
	for line in plasmid_finder:
		start = int(line.split('\t')[6])
		stop = int(line.split('\t')[7])
		scf = line.split('\t')[0]
		mgetype = 'plasmid'
		update_with_scaffold(start, stop, scf, mgetype)

print('plasflow')
with open('/workdir/users/agk85/'+folder+'/plasmids/metagenomes/' + name + '/' + name + '_plasflow.txt') as plasflow:
	header = plasflow.readline()
	for line in plasflow:
		classification = line.split('\t')[5]
		if 'plasmid' in classification:
			start = 0
			stop = float(line.split('\t')[3])
			scf = line.split('\t')[2]
			mgetype = 'plasmid'
			update_with_scaffold(start, stop, scf, mgetype)


print('fullplasmids')
with open('/workdir/users/agk85/'+folder+'/plasmids/metagenomes/' + name + '/' + name + '_fullplasmids_filter.out','r') as fullplasmids:
	for line in fullplasmids:
		flag = 0
		start = int(line.split('\t')[6])
		stop = int(line.split('\t')[7])
		scf = line.split('\t')[0]
		mgetype = 'plasmid'
		update_with_scaffold(start, stop, scf, mgetype)



immedb = {}
immedb_file = '/workdir/refdbs/ImmeDB/MGE_sequences.fasta'
with open(immedb_file) as immedb_seqs:
        for line in immedb_seqs:
                if line[0] == '>':
                        immeid = line.split(' ')[0].split('>')[1]
                        element = line.strip().split(' ')[1]
                        immedb[immeid] = element

print('imme')
with open('/workdir/users/agk85/'+folder+'/imme/metagenomes/'  + name + '/' + name + '_imme_filter.out','r') as imme:
	for line in imme:
		start = int(line.split('\t')[6])
		stop = int(line.split('\t')[7])
		scf = line.split('\t')[0]
		hit = line.split('\t')[1]
		if 'Transposon' in immedb[hit]:
			mgetype = 'transposon'
			update_with_scaffold(start, stop, scf, mgetype)
		if 'Prophage' in immedb[hit]:
			mgetype = 'phage'
			update_with_scaffold(start, stop, scf, mgetype)
		if 'IME' in immedb[hit]:
			mgetype = 'ice'
			update_with_scaffold(start, stop, scf, mgetype)
		if 'ICE' in immedb[hit]:
			mgetype = 'ime'
			update_with_scaffold(start, stop, scf, mgetype)
		else:
			continue


print('plasmid_pfams')
with open('/workdir/users/agk85/'+folder+'/plasmids/metagenomes/' + name + '/' + name + '_plasmid_pfam.txt','r') as plasmidpfams:
	for line in plasmidpfams:
		prot = line.split('\t')[0]
		mgetype = 'plasmid'
		seq = protdict[prot]
		add_mge(mgedict, prot, mgetype, seq)

print('mobile_pfams')
with open('/workdir/users/agk85/'+folder+'/annotation/metagenomes/' +name + '_mobilegenes.txt','r') as mobilepfams:
	for line in mobilepfams:
		prot = line.split('\t')[0]
		mgetype = 'mobile'
		seq = protdict[prot]
		add_mge(mgedict, prot, mgetype, seq)


lh = 0
lt = 0
ht = 0
lht = 0
count = 0
with open('/workdir/users/agk85/'+folder+'/mobile/metagenomes/' + name + '_mge.fna', 'w') as mge_outfile:
	for key,entry in mgedict.iteritems():
		count += 1
		plasmid = '0'
		phage = '0'
		transposon = '0'
		ice = '0'
		ime = '0'
		arg = '0'
		mobile = '0'
		types = set(mgedict[key]['type'])
		seq = mgedict[key]['seq']
		if 'plasmid' in types:
			plasmid = '1'
		if 'phage' in types:
			phage = '1'
		if 'transposon' in types:
			transposon = '1'
		if 'ice' in types:
			ice = '1'
		if 'ime' in types:
			ime = '1'
		if 'arg' in types:
			arg = '1'
		if 'mobile' in types:
			mobile= '1'
		mge_outfile.write('>' + key + ' ' + plasmid + '|' + phage + '|' + transposon +'|' + ice+ '|' + ime + '|' + arg +'|' + mobile + '\n' + seq + '\n')
		if (int(plasmid) + int(phage) + int(transposon)) == 3:
			lht += 1
		if (int(plasmid) + int(phage)) == 2:
			lh += 1
		if (int(plasmid) + int(transposon)) == 2:
			lt += 1
		if (int(phage) + int(transposon)) == 2:
			ht += 1	 

print("lht:",float(lht)/count)
print("lh:",float(lh)/count)
print("lt:",float(lt)/count)
print("ht:",float(ht)/count)

