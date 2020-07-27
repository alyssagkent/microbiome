#$ -S /bin/bash
#$ -N prot_seqs
#$ -e /workdir/users/agk85/CDC2/logs/protseqs_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/protseqs_$JOB_ID.err
#$ -t 1-43
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=2G
#$ -q short.q@cbsubrito2

# Goal is to get scfs from scflist file created by get_card_scfs.sh  
SCRIPT=/workdir/users/agk85/CDC2/scripts/mge_scripts/protseqs_ofinterest7.py
FOLDER=CDC2
WRK=/workdir/users/agk85/${FOLDER}
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`


LOG=/workdir/users/agk85/${FOLDER}/logs/mgeseqs_${NAME}.log
echo start protseqs_ofinterest ${NAME} `date` >> $LOG 2>&1
echo $SCRIPT >> $LOG 2>&1
#run python script
python $SCRIPT $NAME $FOLDER >> $LOG 2>&1

echo end protseqs_ofinterest ${NAME} `date` >> $LOG 2>&1

#after wards

#grep '>' /workdir/users/agk85/CDC/tables/metagenomes3/B*_mge.fna > /workdir/users/agk85/CDC/tables/metagenomes3/gene_ids.txt
