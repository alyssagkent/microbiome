WRK=/workdir/users/agk85/CDC2

#FILES
CONTACT_THRESH=2
GENECLUSTER=$WRK/args/args_99_nr.fna.clstr.tbl
GENETYPE=arg
CONNECTIONS=$WRK/bins/connections_${GENETYPE}_org_all_${CONTACT_THRESH}.txt
SCRIPT=$WRK/scripts/timelapse_scripts/metaphlan_v_gene.py
BINTABLE=$WRK/das/all_bintables.txt
time python $SCRIPT -b $BINTABLE -c $CONNECTIONS -m $CONTACT_THRESH -gc $GENECLUSTER -g $GENETYPE
Rscript $WRK/scripts/timelapse_scripts/plotting_metaphlan_v_gene.R $CONTACT_THRESH $GENETYPE
