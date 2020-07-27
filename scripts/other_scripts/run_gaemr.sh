#$ -S /bin/bash
#$ -N gaemr_CDC2
#$ -V
#$ -o /workdir/users/agk85/CDC2/logs/gaemr_$JOB_ID.out #####
#$ -e /workdir/users/agk85/CDC2/logs/gaemr_$JOB_ID.err #####
#$ -wd /workdir/users/agk85/CDC2/logs #Your working directory
#$ -pe parenv 4
#$ -l h_vmem=100G
#$ -t 23-32 ##change this
#$ -q long.q@cbsubrito2


#Set dirs
FOLDER=CDC2
WRK=/workdir/users/agk85/${FOLDER}
DESIGN_FILE=${WRK}/MetaDesign.txt
NAME=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)

SCF=$WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta
export PATH=/programs/GAEMR-1.0.1/bin:$PATH
export PYTHONPATH=/programs/GAEMR-1.0.1/

echo $SAMPLE `date` Gaemr start
OUT=$WRK/gaemr/metagenomes/${NAME}
if [ ! -d $OUT ]; then mkdir -p $OUT; fi
cd $OUT
/usr/local/bin/python2.7 /programs/GAEMR-1.0.1/bin/GAEMR.py -s $SCF -t 4 -m IDBA_UD -o ${NAME} --force

echo $SAMPLE `date` Gaemr complete

