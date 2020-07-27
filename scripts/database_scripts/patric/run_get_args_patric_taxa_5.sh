WRK=/workdir/users/agk85/CDC2

#FILES
CONTACT_THRESH=5
GENE=arg
GENECLUSTER=$WRK/args/args_99_nr.fna.clstr
SCRIPT=$WRK/scripts/database_scripts/patric/get_genes_onlyhic_taxa.py
SCRIPT2=$WRK/scripts/database_scripts/patric/ARG_comparisons_taxonomies_patric.py
SCRIPT3=$WRK/scripts/database_scripts/patric/ARG_comparisons_taxonomies_patric_support.py

CONNECTIONS=$WRK/bins_hicsupport/connections_${GENE}_org_all_${CONTACT_THRESH}.txt
OUTFILE=$WRK/bins/patric_figures/${GENE}_${CONTACT_THRESH}_base_hic_taxonomies.txt

time python $SCRIPT -c $CONNECTIONS -gc $GENECLUSTER -o $OUTFILE
time python $SCRIPT2 $CONTACT_THRESH
time python $SCRIPT3 $CONTACT_THRESH

#plotting one second each
time Rscript $WRK/scripts/database_scripts/patric/arg_patric_plotting.R $CONTACT_THRESH $GENE onlyhic patric
time Rscript ~/agk/CDC2/scripts/database_scripts/patric/arg_patric_plotting.R $CONTACT_THRESH $GENE support patric
