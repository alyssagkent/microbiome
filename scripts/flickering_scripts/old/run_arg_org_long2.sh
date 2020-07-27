#$ -S /bin/bash
#$ -N mge_long
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/mge_long_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/mge_long_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/scripts
#$ -l h_vmem=100G
#$ -q long.q@cbsubrito2

WRK=/workdir/users/agk85/CDC2

#preparations
#combined das file
#cat $WRK/das/*/*_DASTool_scaffolds2bin.txt > $WRK/das/all_DASTool_scaffolds2bin.txt

#Checkm
#cat $WRK/das/*/checkm_lineage/*.stats > $WRK/das/all.stats

#cvb file
#head -1 $WRK/bins/B314-1_das_1_contigs_v_bins.txt > $WRK/bins/all_das_1_contigs_v_bins_all.txt
#tail -n +2 -q $WRK/bins/*_das_1_contigs_v_bins.txt >> $WRK/bins/all_das_1_contigs_v_bins_all.txt

#kraken
#cat $WRK/das/*/kraken/*_all_kraken.txt > $WRK/das/all_kraken.txt

#FILES
CONTACT_THRESH=5
SCRIPT=$WRK/scripts/arg_scripts/arg_org_long.py
DAS=$WRK/das/all_DASTool_scaffolds2bin.txt
CVB=$WRK/bins/all_das_1_contigs_v_bins_all.txt
CHECKM=$WRK/das/all.stats
KRAKEN_BINS=$WRK/das/all_kraken.txt
GENECLUSTER=$WRK/args/args_99_nr.fna.clstr.tbl
GENE=arg
#GENECLUSTER=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
#GENE=mge
python $SCRIPT -b $DAS -i $CVB -m $CONTACT_THRESH -c $CHECKM -gc $GENECLUSTER -g $GENE -k $KRAKEN_BINS

