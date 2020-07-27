WRK=/workdir/users/agk85/CDC2

#FILES
CONTACT_THRESH=2
GENE=arg
GENECLUSTER=$WRK/args/args_99_nr.fna.clstr.tbl
SCRIPT=$WRK/scripts/diversity_scripts/diversity_v_connectivity.py
SCRIPT2=$WRK/scripts/diversity_scripts/plotting_diversity_connectivity.R

CONNECTIONS=$WRK/bins/connections_${GENE}_org_all_${CONTACT_THRESH}.txt
OUTFILE=$WRK/bins/diversity_figures/${GENE}_${CONTACT_THRESH}_connection_counts.txt
#GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
#GENE=mge
echo get relevant gene taxonomies
time python $SCRIPT -c $CONNECTIONS -gc $GENECLUSTER -o $OUTFILE

SCRIPT2=$WRK/scripts/diversity_scripts/plotting_diversity_connectivity.R
Rscript $SCRIPT2 $CONTACT_THRESH $GENE
SCRIPT3=$WRK/scripts/diversity_scripts/plotting_diversity_connectivity_levels.R
Rscript $SCRIPT3 $CONTACT_THRESH $GENE
