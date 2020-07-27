WRK=/workdir/users/agk85/CDC2

#FILES
CONTACT_THRESH=5
SCRIPT=$WRK/scripts/supportive_scripts/connections_supportive.py

CVB=$WRK/bins_hicsupport/all_das_1_contigs_v_bins_all.txt
BINTABLE=$WRK/das/all_bintables.txt
CHECKM=$WRK/das/all.stats
GENECLUSTER=$WRK/args/args_99_nr.fna.clstr.tbl
GENE=arg
#GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
#GENE=mge
time python $SCRIPT -b $BINTABLE -i $CVB -m $CONTACT_THRESH -c $CHECKM -gc $GENECLUSTER -g $GENE
