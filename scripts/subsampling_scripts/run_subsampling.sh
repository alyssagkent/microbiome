#$ -S /bin/bash
#$ -N contig_vs_bin
#$ -V
#$ -o /workdir/users/agk85/CDC2/logs/contig_vs_bin_$JOB_ID.out #####
#$ -e /workdir/users/agk85/CDC2/logs/contig_vs_bin_$JOB_ID.err #####
#$ -wd /workdir/users/agk85/CDC2/logs #Your working directory
#$ -pe parenv 1
#$ -l h_vmem=20G
#$ -t 1-43  ##change this
#$ -q long.q@cbsubrito2

#Set dirs
FOLDER=CDC2
WRK=/workdir/users/agk85/${FOLDER}
DESIGN_FILE=${WRK}/MetaDesign.txt
NAME=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)

OUT=$WRK/resampling
if [ ! -d $OUT ]; then mkdir -p $OUT; fi
cd $OUT

SCRIPT=$WRK/scripts/subsampling_scripts/subsampling.py


############TRANS############ARG###################
#FILES
DAS=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt
OUTFILE1=${OUT}/${NAME}_subsampling_2_arg_connections.txt
OUTFILE2=${OUT}/${NAME}_subsampling_5_arg_connections.txt
OUTFILE3=${OUT}/${NAME}_subsampling_2_arg_conbin.txt
OUTFILE4=${OUT}/${NAME}_subsampling_5_arg_conbin.txt
WINDOW=500
HIC=$WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_trans_primary_0_noeuks.txt
COMBO_TABLE=$WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
CHECKM=$WRK/das/${NAME}/checkm_lineage/${NAME}.stats
ARGS=$WRK/args/args_99_nr.fna.clstr.tbl
#MGE=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
KRAKEN_BINS=$WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt
time python $SCRIPT -b $DAS -l $HIC -c $CHECKM -gc $ARGS -k $KRAKEN_BINS -o1 $OUTFILE1 -o2 $OUTFILE2 -o3 $OUTFILE3 -o4 $OUTFILE4 -w $WINDOW

############TRANS############MGE###################
#FILES
DAS=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt
OUTFILE1=${OUT}/${NAME}_subsampling_2_mge_connections.txt
OUTFILE2=${OUT}/${NAME}_subsampling_5_mge_connections.txt
OUTFILE3=${OUT}/${NAME}_subsampling_2_mge_conbin.txt
OUTFILE4=${OUT}/${NAME}_subsampling_5_mge_conbin.txt
WINDOW=500
HIC=$WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_trans_primary_0_noeuks.txt
COMBO_TABLE=$WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
CHECKM=$WRK/das/${NAME}/checkm_lineage/${NAME}.stats
#ARGS=$WRK/args/args_99_nr.fna.clstr.tbl
MGE=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
KRAKEN_BINS=$WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt
time python $SCRIPT -b $DAS -l $HIC -c $CHECKM -gc $MGE -k $KRAKEN_BINS -o1 $OUTFILE1 -o2 $OUTFILE2 -o3 $OUTFILE3 -o4 $OUTFILE4 -w $WINDOW

############ALLVP############ARG###################
#FILES
DAS=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt
OUTFILE1=${OUT}/${NAME}_subsampling_2_arg_connections_vp.txt
OUTFILE2=${OUT}/${NAME}_subsampling_5_arg_connections_vp.txt
OUTFILE3=${OUT}/${NAME}_subsampling_2_arg_conbin_vp.txt
OUTFILE4=${OUT}/${NAME}_subsampling_5_arg_conbin_vp.txt
WINDOW=50000
HIC=$WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_allValidPairs
COMBO_TABLE=$WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
CHECKM=$WRK/das/${NAME}/checkm_lineage/${NAME}.stats
ARGS=$WRK/args/args_99_nr.fna.clstr.tbl
#MGE=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
KRAKEN_BINS=$WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt
time python $SCRIPT -b $DAS -l $HIC -c $CHECKM -gc $ARGS -k $KRAKEN_BINS -o1 $OUTFILE1 -o2 $OUTFILE2 -o3 $OUTFILE3 -o4 $OUTFILE4 -w $WINDOW

############ALLVP############MGE###################
#FILES
DAS=$WRK/das/${NAME}/${NAME}_DASTool_scaffolds2bin.txt
OUTFILE1=${OUT}/${NAME}_subsampling_2_mge_connections_vp.txt
OUTFILE2=${OUT}/${NAME}_subsampling_5_mge_connections_vp.txt
OUTFILE3=${OUT}/${NAME}_subsampling_2_mge_conbin_vp.txt
OUTFILE4=${OUT}/${NAME}_subsampling_5_mge_conbin_vp.txt
WINDOW=50000
HIC=$WRK/hicpro/output/${NAME}_output/hic_results/data/${NAME}/${NAME}_allValidPairs
COMBO_TABLE=$WRK/combo_tables/metagenomes/${NAME}_master_scf_table.txt
CHECKM=$WRK/das/${NAME}/checkm_lineage/${NAME}.stats
#ARGS=$WRK/args/args_99_nr.fna.clstr.tbl
MGE=$WRK/mobile/metagenomes/mge_99_nr.fna.clstr.tbl
KRAKEN_BINS=$WRK/das/${NAME}/kraken/${NAME}_all_kraken_weighted.txt
time python $SCRIPT -b $DAS -l $HIC -c $CHECKM -gc $MGE -k $KRAKEN_BINS -o1 $OUTFILE1 -o2 $OUTFILE2 -o3 $OUTFILE3 -o4 $OUTFILE4 -w $WINDOW

