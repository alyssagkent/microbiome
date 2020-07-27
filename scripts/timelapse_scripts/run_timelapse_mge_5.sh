WRK=/workdir/users/agk85/CDC2

#FILES
CONTACT_THRESH=5
GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
GENE=mge
CONNECTIONS=$WRK/bins/connections_${GENE}_org_all_${CONTACT_THRESH}.txt
SCRIPT=$WRK/scripts/timelapse_scripts/timelapse.py
CVB=$WRK/bins/all_das_1_contigs_v_bins_all.txt
BINTABLE=$WRK/das/all_bintables.txt
time python $SCRIPT -b $BINTABLE -i $CVB -c $CONNECTIONS -m $CONTACT_THRESH -gc $GENECLUSTER -g $GENE

Rscript $WRK/scripts/timelapse_scripts/plotting_timelapse.R $CONTACT_THRESH $GENETYPE
