SCRIPT=/workdir/users/agk85/CDC2/scripts/org_org_scripts/taxa_taxa_prep_healthy_sick.py
SCRIPT2=/workdir/users/agk85/CDC2/scripts/org_org_scripts/plotting_taxa_taxa_taxonpresence2_healthy_sick.R
GENE=arg
THRESH=2
python ${SCRIPT} -i /workdir/users/agk85/CDC2/bins/histogram_figures/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_healthysick.txt

GENE=arg
THRESH=5
python ${SCRIPT} -i /workdir/users/agk85/CDC2/bins/histogram_figures/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_healthysick.txt

GENE=mge
THRESH=2
python ${SCRIPT} -i /workdir/users/agk85/CDC2/bins/histogram_figures/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_healthysick.txt

GENE=mge
THRESH=5
python ${SCRIPT} -i /workdir/users/agk85/CDC2/bins/histogram_figures/Together_das_${THRESH}_${GENE}taxa_together.txt -o /workdir/users/agk85/CDC2/bins/${GENE}_das_${THRESH}_taxa_taxa_healthysick.txt


SCRIPT2=/workdir/users/agk85/CDC2/scripts/org_org_scripts/plotting_taxa_taxa_taxonpresence2_healthy_sick.R
Rscript ${SCRIPT2} /workdir/users/agk85/CDC2/bins/arg_das_2_taxa_taxa_healthysick.txt 2 arg 10
Rscript ${SCRIPT2} /workdir/users/agk85/CDC2/bins/arg_das_5_taxa_taxa_healthysick.txt 5 arg 10
Rscript ${SCRIPT2} /workdir/users/agk85/CDC2/bins/mge_das_2_taxa_taxa_healthysick.txt 2 mge 400
Rscript ${SCRIPT2} /workdir/users/agk85/CDC2/bins/mge_das_5_taxa_taxa_healthysick.txt 5 mge 400

cp ~/agk/CDC2/bins/org_org_figures/Org*hs.pdf ~/CDC_figures/org_org
