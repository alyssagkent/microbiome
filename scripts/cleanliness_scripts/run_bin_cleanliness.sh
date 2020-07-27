#$ -S /bin/bash
#$ -N bin_cleanliness
#$ -V
#$ -o /workdir/users/agk85/CDC2/logs/bin_cleanliness_$JOB_ID.out #####
#$ -e /workdir/users/agk85/CDC2/logs/bin_cleanliness_$JOB_ID.err #####
#$ -wd /workdir/users/agk85/CDC2/logs #Your working directory
#$ -pe parenv 1
#$ -l h_vmem=20G
#$ -t 1-43  ##change this
#$ -q short.q@cbsubrito2

#Set dirs
FOLDER=CDC2
WRK=/workdir/users/agk85/${FOLDER}
DESIGN_FILE=${WRK}/MetaDesign.txt
NAME=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)

OUT=$WRK/bins
if [ ! -d $OUT ]; then mkdir -p $OUT; fi
cd $OUT

#FILES
SCRIPT=$WRK/scripts/cleanliness_scripts/bin_cleanliness.py
DAS=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt
OUTFILE=$WRK/bins/${NAME}_bin_cleanliness.txt
HIC=$WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_allValidPairs
COMBO_TABLE=$WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
CHECKM=$WRK/das/${NAME}/checkm_lineage/${NAME}.stats
KRAKEN_BINS=$WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt
time python $SCRIPT -b $DAS -o $OUTFILE -l $HIC -t $COMBO_TABLE -c $CHECKM -k $KRAKEN_BINS

