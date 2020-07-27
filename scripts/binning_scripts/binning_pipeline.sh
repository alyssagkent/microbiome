


#run maxbin
qsub_run_maxbin.sh
#run metabat
#first you have to map the reads to the scaffolds
qsub snake_bwa_scaffold_qsub.sh 

screen
bash run_metabat_manual.sh
#run Concoct
#this maps the patients reads to each of its timepoints scaffolds (one timepoint at a time)
qsub snake_bwa_scaffold_t1_qsub.sh
qsub snake_bwa_scaffold_t2_qsub.sh
qsub snake_bwa_scaffold_t3_qsub.sh
qsub snake_bwa_scaffold_t4_qsub.sh
qsub snake_bwa_scaffold_t5_qsub.sh
qsub snake_bwa_scaffold_t6_qsub.sh
qsub snake_bwa_scaffold_t7_qsub.sh


qsub run_concoct_t1.sh
qsub run_concoct_t2.sh
qsub run_concoct_t3.sh
qsub run_concoct_t4.sh
qsub run_concoct_t5.sh
qsub run_concoct_t6.sh
qsub run_concoct_t7.sh

#make sure it has the right range
qsub run_rename_concoct.sh

#run DAS

#run checkm on each of them
qsub run_checkm_maxbin.sh
qsub run_checkm_metabat.sh
qsub run_checkm_concoct_timepoints.sh

#this one handles the stats well..
qsub run_checkm_das.sh

#combine all of the stats files and put them into the maxbin folder (convenience)

#maxbin stats

#metabat stats 

#concoct stats
cd /workdir/users/agk85/CDC2/concoct/
for file in */checkm_lineage/*.stats; do  echo $file; sample=$(echo `basename "$file"` | cut -d'.' -f1); echo $sample; awk -v samp="${sample}_" '{print samp $0}' $file> ${file}.named; done
cat */checkm_lineage/*.named > ~/agk/CDC2/maxbin/concoct_checkm_stats.txt

#das stats
cat ~/agk/CDC2/das/*/checkm_lineage/*.stats > ~/agk/CDC2/maxbin/das_withconcoct_checkm_stats.txt

Rscript binner_comparison_checkm.R


cd /workdir/users/agk85/CDC2/das/
cat */checkm_lineage/*.stats | sort > all_das_checkm.txt
cat */kraken/*.contamination | sort > all_kraken_besttaxid.txt
join -j 1 all_das_checkm.txt all_kraken_besttaxid.txt -t $'\t' > all_das_checkm_kraken.txt
cp all_das_checkm_kraken.txt ../bins/

