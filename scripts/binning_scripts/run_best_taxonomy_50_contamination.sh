#$ -S /bin/bash
#$ -N kraken_taxonomy
#$ -e /workdir/users/agk85/CDC2/logs/kraken_taxonomy_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/kraken_taxonomy_$JOB_ID.out
#$ -t 1-36
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito2

WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/concoct/MetaDesign_concoct.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`


for file in $WRK/das/${NAME}/kraken/*report.txt ;
do 
echo $file;
outfile=${file}.besttaxid;
python /workdir/users/agk85/CDC2/scripts/best_taxonomy_50_contamination.py -i $file -o $outfile;
done


