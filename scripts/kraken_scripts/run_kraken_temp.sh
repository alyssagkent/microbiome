#$ -S /bin/bash
#$ -N kraken_bins
#$ -e /workdir/users/agk85/CDC2/logs/kraken_bins_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/kraken_bins_$JOB_ID.out
#$ -t 1-32
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -l h_vmem=1G
#$ -q short.q@cbsubrito2

# Goal is to bin contigs into single species from idba generated scaffolds

WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

OUT=$WRK/das/${NAME}/kraken/
if [ ! -d $OUT ]; then mkdir -p $OUT; fi


#merge all of the besttaxid files together with the name
for file in $OUT/*.besttaxid; do name=`basename "$file"`; sed -i "s/^/$name\t/" $file; done
cat $OUT/*besttaxid > $OUT/${NAME}_all_kraken.txt


