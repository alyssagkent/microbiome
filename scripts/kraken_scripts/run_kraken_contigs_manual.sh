WRK=/workdir/users/agk85/CDC2
for SGE_TASK_ID in {1..16};
do
#list of all the bins in the DAS folders:
#ls */*DASTool_bins/*.fa | cut -d'/' -f3 > bins.txt
DESIGN_FILE=$WRK/MetaDesign_newmgms.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

#make the bin file

mkdir $WRK/kraken/${NAME}

echo start kraken for $NAME
/workdir/users/agk85/tools/krakenuniq/krakenuniq --db /workdir/users/agk85/tools/krakenuniq/DB_subset --threads 32 --fasta-input $WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta --report-file $WRK/kraken/${NAME}/${NAME}.report.txt > ${WRK}/kraken/${NAME}/${NAME}.kraken.txt 

/workdir/users/agk85/tools/krakenuniq/krakenuniq-translate --db /workdir/users/agk85/tools/krakenuniq/DB_subset --mpa-format ${WRK}/kraken/${NAME}/${NAME}.kraken.txt  > ${WRK}/kraken/${NAME}/${NAME}.kraken.taxonomy.txt

echo end kraken for $NAME
done
