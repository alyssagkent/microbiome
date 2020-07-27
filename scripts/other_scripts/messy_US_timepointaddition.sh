
for file in *; do echo $file; sed -i 's/US2/US2-1/g' $file; done

mv US2.gff3 US2-1.gff3
US2_excise.fasta US2-1_excise.fasta
US2_phagefinder.txt US2-1_phagefinder.txt 
US2_proteins.faa 
US2_proteins.fna
US2_scaffold.fasta
US2_scaffold.fasta.fai
US2_scfids.txt
