
LOG=/workdir/users/agk85/CDC/tables/log/cdhit_log.txt

PERCENT=0.95
NUM=95
/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i mobile_genes.fna -o mge_${NUM}_nr.fna -n 6 -c $PERCENT -M 30000 -T 15 -d 0 -s .9 >> $LOG 2>&1



