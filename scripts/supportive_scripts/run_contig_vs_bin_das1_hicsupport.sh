#$ -S /bin/bash
#$ -N contig_vs_bin
#$ -V
#$ -o /workdir/users/agk85/CDC2/logs/contig_vs_bin_$JOB_ID.out #####
#$ -e /workdir/users/agk85/CDC2/logs/contig_vs_bin_$JOB_ID.err #####
#$ -wd /workdir/users/agk85/CDC2/logs #Your working directory
#$ -l h_vmem=10G
#$ -t 1-43  ##change this
#$ -q short.q@cbsubrito2

#working through 22-43

#Set dirs
FOLDER=CDC2
CONTACT_THRESH=1
WRK=/workdir/users/agk85/${FOLDER}
DESIGN_FILE=${WRK}/MetaDesign.txt
NAME=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)

SCF=$WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta

echo $SAMPLE `date` CvB start
OUT=$WRK/bins_hicsupport
if [ ! -d $OUT ]; then mkdir -p $OUT; fi
cd $OUT

#FILES
SCRIPT=$WRK/scripts/supportive_scripts/contig_vs_bin_hicsupport.py
DAS=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt
OUTFILE=${OUT}/${NAME}_das_${CONTACT_THRESH}_contigs_v_bins.txt
HIC=$WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_trans_primary_0_ncol_withexcise_noeuks_normalize_${CONTACT_THRESH}.txt
COMBO_TABLE=$WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
CHECKM=$WRK/das/${NAME}/checkm_lineage/${NAME}.stats
ARGS=$WRK/args/args_99_nr.fna.clstr.tbl
MGE=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
KRAKEN_BINS=$WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt
KRAKEN_CONTIGS=$WRK/kraken/${NAME}/${NAME}.kraken.taxonomy.txt
time python $SCRIPT -b $DAS -o $OUTFILE -l $HIC -t $COMBO_TABLE -m $CONTACT_THRESH -c $CHECKM -ar $ARGS -mge $MGE -k $KRAKEN_BINS -kc $KRAKEN_CONTIGS

#python $WRK/scripts/contig_histogram_prep.py  -o $WRK/bins/${NAME}_contig_taxa_network_${CONTACT_THRESH}_fgs  -i $WRK/bins/${NAME}_network_${CONTACT_THRESH}_contigs_v_bins.txt  -c $WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
echo $SAMPLE `date` CvB complete
