#$ -S /bin/bash
#$ -N gaemr_gettaxa
#$ -V
#$ -t 19-32
#$ -e /workdir/users/agk85/CDC2/logs/gaemr_gettaxa_$JOB_ID.e
#$ -o /workdir/users/agk85/CDC2/logs/gaemr_gettaxa_$JOB_ID.o
#$ -wd /workdir/users/agk85/CDC2
#$ -l h_vmem=2G
#$ -pe parenv 1
#$ -q long.q@cbsubrito2

FOLDER=CDC2
WRK=/workdir/users/agk85/${FOLDER}
DESIGN_FILE=$WRK/MetaDesign.txt #changed from rerun
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

echo combo_table ${NAME} `date`
python  ~/agk/US3/scripts/gaemr_gettaxa_percent.py -i /workdir/users/agk85/${FOLDER}/gaemr/metagenomes/ -n $NAME -s /workdir/users/agk85/${FOLDER}/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta
