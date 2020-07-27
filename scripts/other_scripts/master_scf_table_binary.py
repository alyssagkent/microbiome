#convert master_scf_table to binary
import sys
folder = sys.argv[1]
name = sys.argv[2]

#name = 'B320-5'
with open("/workdir/users/agk85/" + folder + "/combo_tables/metagenomes/" + name + '_master_scf_table.txt','r') as infile:
	with open("/workdir/users/agk85/" + folder + "/combo_tables/metagenomes/" + name + '_master_scf_table_binary.txt','w') as outfile:
		header = infile.readline()
		#outfile.write(header)
		headerfields =header.strip().split('\t') 
		outfile.write(headerfields.pop(0))
		for item in headerfields:
			outfile.write('\t' + item.split('(')[0])
		outfile.write('\n')
		for line in infile:
			fields = line.strip().split('\t')
			id = fields.pop(0)
			outfile.write(id)
			for field in fields:
				if field != '.':
					outfile.write('\t1')
				else:
					outfile.write('\t0')
			outfile.write('\n')
