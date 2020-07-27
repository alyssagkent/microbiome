#$ -S /bin/bash
#$ -N remove_mobile
#$ -e /workdir/users/agk85/CDC2/logs/remove_mobile_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/remove_mobile_$JOB_ID.out
#$ -t 1-37
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito2

WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

OUT=${WRK}/das/${NAME}/${NAME}_DASTool_non_mobile_bins/
if [ ! -d $OUT ]; then mkdir -p $OUT; fi


for file in $WRK/das/${NAME}/${NAME}_DASTool_bins/*.fa;
do 
echo $file;
outfile="${file/\.fa/_non_mobile.fa}"
outfile2="${outfile/DASTool_bins/DASTool_non_mobile_bins}"
python /workdir/users/agk85/CDC2/scripts/non_mobile_clusters.py -i $file -o $outfile2 -c ${WRK}/combo_tables;

done


