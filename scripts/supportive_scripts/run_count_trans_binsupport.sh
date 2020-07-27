#$ -S /bin/bash
#$ -N hic_support_trans
#$ -V
#$ -o /workdir/users/agk85/CDC2/logs/hic_support_trans_$JOB_ID.out #####
#$ -e /workdir/users/agk85/CDC2/logs/hic_support_trans_$JOB_ID.err #####
#$ -wd /workdir/users/agk85/CDC2/logs #Your working directory
#$ -l h_vmem=1G
#$ -t 1-43  ##change this
#$ -q short.q@cbsubrito2

#working through 22-43

#Set dirs
FOLDER=CDC2
WRK=/workdir/users/agk85/${FOLDER}
DESIGN_FILE=${WRK}/MetaDesign.txt
NAME=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)

OUT=$WRK/bins_hicsupport
if [ ! -d $OUT ]; then mkdir -p $OUT; fi
cd $OUT

#FILES
SCRIPT=$WRK/scripts/supportive_scripts/count_trans_binsupport.py
DAS=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt
OUTFILE=${OUT}/${NAME}_trans_binsupport.txt
HIC=$WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_allValidPairs
COMBO_TABLE=$WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
time python $SCRIPT -b $DAS -o $OUTFILE -l $HIC -t $COMBO_TABLE 

