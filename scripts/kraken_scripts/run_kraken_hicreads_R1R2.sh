#$ -S /bin/bash
#$ -N kraken_reads
#$ -e /workdir/users/agk85/CDC2/logs/kraken_reads_sep_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/kraken_reads_sep_$JOB_ID.out
#$ -t 20-22
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=300G
#$ -q long.q@cbsubrito2

# Goal is to bin contigs into single species from idba generated scaffolds

WRK=/workdir/users/agk85/CDC2

#list of all the bins in the DAS folders:
#ls */*DASTool_bins/*.fa | cut -d'/' -f3 > bins.txt
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

#make the bin file

OUT=$WRK/kraken/${NAME}
if [ ! -d $OUT ]; then mkdir -p $OUT; fi

echo start kraken for $NAME

/workdir/users/agk85/tools/krakenuniq/krakenuniq --db /workdir/users/agk85/tools/krakenuniq/DB_subset --threads 32 --fastq-input --report-file $WRK/kraken/${NAME}/${NAME}.hicreads.R1.report.txt /workdir/data/CDC/hic/merged/${NAME}hic.1.fastq > ${WRK}/kraken/${NAME}/${NAME}.hicreads.R1.kraken.txt

/workdir/users/agk85/tools/krakenuniq/krakenuniq --db /workdir/users/agk85/tools/krakenuniq/DB_subset --threads 32 --fastq-input --report-file $WRK/kraken/${NAME}/${NAME}.reads.R2.report.txt /workdir/data/CDC/hic/merged/${NAME}hic.2.fastq > ${WRK}/kraken/${NAME}/${NAME}.hicreads.R2.kraken.txt

/workdir/users/agk85/tools/krakenuniq/krakenuniq-translate --db /workdir/users/agk85/tools/krakenuniq/DB_subset --mpa-format ${WRK}/kraken/${NAME}/${NAME}.hicreads.R1.kraken.txt  > ${WRK}/kraken/${NAME}/${NAME}.hicreads.R1.kraken.taxonomy.txt
/workdir/users/agk85/tools/krakenuniq/krakenuniq-translate --db /workdir/users/agk85/tools/krakenuniq/DB_subset --mpa-format ${WRK}/kraken/${NAME}/${NAME}.hicreads.R2.kraken.txt  > ${WRK}/kraken/${NAME}/${NAME}.hicreads.R2.kraken.taxonomy.txt

echo end kraken for $NAME
