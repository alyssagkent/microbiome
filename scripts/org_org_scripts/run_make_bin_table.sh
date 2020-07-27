for SGE_TASK_ID in {1..43};
do
WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

SCRIPT=$WRK/scripts/org_org_scripts/make_bin_table.py
CHECKM=$WRK/das/${NAME}/checkm_lineage/${NAME}.stats
KRAKEN_BIN=$WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt
OUTFILE=$WRK/das/${NAME}/${NAME}_bintable.txt
python $SCRIPT -c $CHECKM -k $KRAKEN_BIN -o $OUTFILE
done
