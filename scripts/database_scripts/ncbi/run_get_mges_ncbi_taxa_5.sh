#!/bin/bash

WRK=/workdir/users/agk85/CDC2

#FILES
CONTACT_THRESH=5
GENE=mge
GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr
SCRIPT=$WRK/scripts/database_scripts/ncbi/get_machinery_genes_onlyhic_taxa.py
SCRIPT2=$WRK/scripts/database_scripts/ncbi/phage_machinery_comparisons_taxonomies_ncbi.py
SCRIPT3=$WRK/scripts/database_scripts/ncbi/phage_machinery_comparisons_taxonomies_ncbi_support.py

CONNECTIONS=$WRK/bins_hicsupport/connections_${GENE}_org_all_${CONTACT_THRESH}.txt
OUTFILE=$WRK/bins/patric_figures/${GENE}_${CONTACT_THRESH}_base_hic_taxonomies.txt
#GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
#GENE=mge
time python $SCRIPT -c $CONNECTIONS -gc $GENECLUSTER -o $OUTFILE -m $CONTACT_THRESH
time python $SCRIPT2 $CONTACT_THRESH
time python $SCRIPT3 $CONTACT_THRESH

#plotting one second each
time Rscript ${WRK}/scripts/database_scripts/ncbi/phage_machinery_plotting.R $CONTACT_THRESH $GENE onlyhic ncbi
time Rscript ${WRK}/scripts/database_scripts/ncbi/phage_machinery_plotting.R $CONTACT_THRESH $GENE support ncbi
