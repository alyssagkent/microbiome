#counts number of phage that are relinked from hic trans data
#this should be run in /workdir/users/agk85/CDC2/hicpro/outputs/
phage = []
joint = 0
joint2 = 0
joint5 = 0
with open('phage_hic.txt') as infile:
	for line in infile:
		c1 = line.split(':')[1].split('\t')[0]
		c2 = line.split('\t')[1]
		count =float(line.strip().split('\t')[2])
		if 'phage' in c1:
			phage.append(c1)
			simplified = c1.split('_')[0] + '_scaffold_' + c1.split('_')[1].split('|')[0]
			if c2 == simplified:
				joint+=1
				if count>1:
					joint2+=1
				if count>4:
					joint5+=1
		if 'phage' in c2:
			phage.append(c2)
			simplified = c2.split('_')[0] + '_scaffold_' + c2.split('_')[1].split('|')[0]
			if c1 == simplified:
				joint+=1
				if count>1:
					joint2+=1
				if count>4:
					joint5+=1

print(len(set(phage)))
print(joint)
print(joint2)
print(joint5)
