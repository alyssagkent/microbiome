WRK=/workdir/users/agk85/CDC2

#FILES
CONTACT_THRESH=2
SCRIPT=$WRK/scripts/arg_scripts/mge_org_long.py
DAS=$WRK/das/all_DASTool_scaffolds2bin.txt
CVB=$WRK/bins/all_das_1_contigs_v_bins_all.txt
CHECKM=$WRK/das/all.stats
KRAKEN_BINS=$WRK/das/all_kraken.txt
#GENECLUSTER=$WRK/args/args_99_nr.fna.clstr.tbl
#GENE=arg
GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
GENE=mge
python $SCRIPT -b $DAS -i $CVB -m $CONTACT_THRESH -c $CHECKM -gc $GENECLUSTER -g $GENE -k $KRAKEN_BINS

