#$ -S /bin/bash
#$ -N mash
#$ -e /workdir/users/agk85/CDC2/logs/mash_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/mash_$JOB_ID.out
#$ -t 1
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=400G
#$ -q long.q@cbsubrito2

#conda activate forDAS
#each of the folders
WRK=/workdir/users/agk85/CDC2
OUT=$WRK/bins/all_bins
#We copied all of the bins here
#cp  $WRK/das/*/*DASTool_bins/*.fa $OUT


#
#cd $OUT
#/programs/mash/mash sketch -s 1000000 -o $OUT *.fa
#/programs/mash/mash dist -t $WRK/bins/all_bins.msh $WRK/bins/all_bins.msh > all_bins_dist.txt

OUT=$WRK/bins/hq_bins
cd $OUT
#/programs/mash/mash sketch -s 1000000 -o ${OUT}_mil *.fa
/programs/mash/mash dist -t ${OUT}_mil.msh ${OUT}_mil.msh > ${OUT}_mil_dist.txt

OUT=$WRK/bins/mq_bins
cd $OUT
#/programs/mash/mash sketch -s 1000000 -o ${OUT}_mil *.fa
/programs/mash/mash dist -t ${OUT}_mil.msh ${OUT}_mil.msh > ${OUT}_mil_dist.txt

OUT=$WRK/bins/quality_bins
cd $OUT
#/programs/mash/mash sketch -s 1000000 -o ${OUT}_mil *.fa
/programs/mash/mash dist -t ${OUT}_mil.msh ${OUT}_mil.msh > ${OUT}_mil_dist.txt
