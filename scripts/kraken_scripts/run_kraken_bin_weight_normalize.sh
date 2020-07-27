#$ -S /bin/bash
#$ -N kraken_bins
#$ -e /workdir/users/agk85/CDC2/logs/kraken_bins_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/kraken_bins_$JOB_ID.out
#$ -t 1-43
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=1G
#$ -q short.q@cbsubrito2

# Goal is to bin contigs into single species from idba generated scaffolds

WRK=/workdir/users/agk85/CDC2

#list of all the bins in the DAS folders:
#ls */*DASTool_bins/*.fa | cut -d'/' -f3 > bins.txt
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

#make the bin file
#ls $WRK/das/${NAME}/${NAME}_DASTool_bins/*.fa | cut -d '/' -f9 > $WRK/das/${NAME}/${NAME}_bins.txt 

OUT=$WRK/das/${NAME}/kraken/
if [ ! -d $OUT ]; then mkdir -p $OUT; fi

#/workdir/users/agk85/tools/krakenuniq/krakenuniq --db /workdir/users/agk85/tools/krakenuniq/DB_subset --threads 10 --preload
file_lines=`cat $WRK/das/${NAME}/${NAME}_bins.txt`
for line in $file_lines; 
do
echo $line
#echo start kraken for DAS $line
#/workdir/users/agk85/tools/krakenuniq/krakenuniq --db /workdir/users/agk85/tools/krakenuniq/DB_subset --threads 4 --fasta-input $WRK/das/${NAME}/${NAME}_DASTool_bins/${line} --report-file $WRK/das/${NAME}/kraken/${line}.report.txt > ${WRK}/das/${NAME}/kraken/${line}.kraken.txt 

#/workdir/users/agk85/tools/krakenuniq/krakenuniq-translate --db /workdir/users/agk85/tools/krakenuniq/DB_subset --mpa-format ${WRK}/das/${NAME}/kraken/${line}.kraken.txt  > ${WRK}/das/${NAME}/kraken/${line}.kraken.taxonomy.txt

BIN=${line%".contigs.fa"}
echo $BIN
#get the best taxonomy for each bin
SCRIPT=$WRK/scripts/kraken_scripts/kraken_bin_weight_normalize.py
INFILE=$WRK/das/${NAME}/kraken/${BIN}.contigs.fa.kraken.taxonomy.txt
FASTA=$WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta
OUTFILE=$WRK/das/${NAME}/kraken/${BIN}.contigs.fa.kraken.weighted_besttaxid.txt
DASFILE=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt

python $SCRIPT -i $INFILE -f $FASTA -o $OUTFILE -d $DASFILE -b $BIN
#cat the besttaxonomy together with the name
echo end kraken for DAS $line
done

#merge all of the besttaxid files together with the name
cat $OUT/*.weighted_besttaxid.txt > $OUT/${NAME}_all_kraken_weighted.txt

