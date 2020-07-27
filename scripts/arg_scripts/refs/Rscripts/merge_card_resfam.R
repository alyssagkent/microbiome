infile1 = read.csv('card_prot_id_type.txt', sep="\t",header=F)
infile2 = read.csv('resfams_prot_id_type.txt', sep="\t",header=F)
merged = merge(infile1, infile2, by="V1",all=T)
write.table(merged, 'arg_prot_id_type.txt')

#aftewards
#sed -i -e 's/"//g' arg_prot_id_type.txt
