
#this grabs all of the mobile genes from all of the different sources including contigs--wraps over the whole contig and grabs all the genes
qsub run_protseqs_ofinterest.sh
#cat all of the fnas into mobile_genes.fn

cat *_mge.fna > mobile_proteins.fna

#cd-hit
bash ~/agk/CDC2/scripts/mge_scripts/run_cdhit_all.sh

#get the cd-hit table
python ~/agk/CDC2/scripts/mge_scripts/cdhit_output_to_table.py CDC2

#get all of the gene_ids
grep '>' ~/agk/CDC2/mobile/metagenomes/mobile_proteins.fna > ~/agk/CDC2/mobile/metagenomes/gene_ids.txt
python ~/agk/CDC2/scripts/mge_scripts/grab_mge_info.py CDC2
sed -i -e "s/'/ prime/g" mge_95_nr.fna.clstr.desc  to make sure the primes are not disrupting R


#prepare machinery genes
cat ~/agk/CDC2/hicpro/output/*_output/hic_results/data/*/*_trans_primary_0_ncol_withexcise_noeuks.txt > ~/agk/CDC2/mobile/metagenomes/CDC+healthy_trans_primary_ncol_0_withexcise_noeuks.txt
