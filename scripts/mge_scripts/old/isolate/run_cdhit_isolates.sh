LOG=/workdir/users/agk85/CDC/tables/log/cdhit_log.txt

PERCENT=0.80
NUM=80

/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i all_genome_proteins.fna -o isolate_${NUM}_nr.fna -n 10 -c $PERCENT -M 30000 -T 15 -d 0 -s .9 >> $LOG 2>&1
#/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i mobile_isolate_genes.fna -o mge_isolate_${NUM}_nr.fna -n 10 -c $PERCENT -M 30000 -T 15 -d 0 -s .9 >> $LOG 2>&1
