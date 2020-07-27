#! /bin/bash

WRK=/workdir/users/agk85/CDC2
GENETYPE=arg
CONTACT_THRESH=5
GENECLUSTER=$WRK/args/args_99_nr.fna.clstr
#GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr

SCRIPT=$WRK/scripts/hgt_scripts/hgt_comparisons_allorg.py
CONNECTIONS=$WRK/bins/connections_${GENETYPE}_org_all_${CONTACT_THRESH}.txt
OUTHANDLE=$WRK/bins/hgt_comparisons/hgt_comparisons_alllevels_allorgs_${GENETYPE}_${CONTACT_THRESH}.txt
OUT=$WRK/bins/hgt_comparisons
if [ ! -d $OUT ]; then mkdir -p $OUT; fi

python $SCRIPT -c $CONNECTIONS -gc $GENECLUSTER -o $OUTHANDLE

Rscript $WRK/scripts/hgt_scripts/hgt_comparisons_plotting.R $CONTACT_THRESH $GENETYPE
