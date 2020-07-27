LOG=/workdir/users/agk85/CDC2/logs/cdhit_log.txt

PERCENT=0.99
NUM=99
/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i /workdir/users/agk85/CDC2/mobile/metagenomes/mobile_proteins.fna -o /workdir/users/agk85/CDC2/mobile/metagenomes/mge_${NUM}_nr.fna -n 10 -c $PERCENT -M 30000 -T 32 -d 0 -s .9 >> $LOG 2>&1

#/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i /workdir/users/agk85/CDC2/mobile/metagenomes/all_proteins.fna -o /workdir/users/agk85/CDC2/mobile/metagenomes/all_${NUM}_nr.fna -n 10 -c $PERCENT -M 30000 -T 32 -d 0 -s .9 >> $LOG 2>&1
