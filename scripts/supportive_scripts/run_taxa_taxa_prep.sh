SCRIPTFOLD=/workdir/users/agk85/CDC2/scripts/org_org_scripts

GENE=arg
THRESH=2
python ${SCRIPTFOLD}/taxa_taxa_prep.py -i /workdir/users/agk85/CDC2/bins/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_counts.txt

GENE=arg
THRESH=5
python ${SCRIPTFOLD}/taxa_taxa_prep.py -i /workdir/users/agk85/CDC2/bins/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_counts.txt

GENE=mge
THRESH=2
python ${SCRIPTFOLD}/taxa_taxa_prep.py -i /workdir/users/agk85/CDC2/bins/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_counts.txt

GENE=mge
THRESH=5
python ${SCRIPTFOLD}/taxa_taxa_prep.py -i /workdir/users/agk85/CDC2/bins/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_counts.txt



Rscript ${SCRIPTFOLD}/plotting_taxa_taxa_taxonpresence2.R /workdir/users/agk85/CDC2/bins/arg_das_2_taxa_taxa_counts.txt 2 arg 5
Rscript ${SCRIPTFOLD}/plotting_taxa_taxa_taxonpresence2.R /workdir/users/agk85/CDC2/bins/arg_das_5_taxa_taxa_counts.txt 5 arg 5
Rscript ${SCRIPTFOLD}/plotting_taxa_taxa_taxonpresence2.R /workdir/users/agk85/CDC2/bins/mge_das_2_taxa_taxa_counts.txt 2 mge 1000
Rscript ${SCRIPTFOLD}/plotting_taxa_taxa_taxonpresence2.R /workdir/users/agk85/CDC2/bins/mge_das_5_taxa_taxa_counts.txt 5 mge 1000
