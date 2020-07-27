#$ -S /bin/bash
#$ -N combo_table_clusters
#$ -e /workdir/users/agk85/CDC2/logs/combo_table_clusters_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/combo_table_clusters_$JOB_ID.out
#$ -t 1-32
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=1G
#$ -q long.q@cbsubrito2

WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

#input files
#combo_tables
INFILE=${WRK}/combo_tables/metagenomes/${NAME}_master_scf_table.txt
OUTFILE=${WRK}/combo_tables/metagenomes/${NAME}_master_scf_cluster_table.txt
DAS=${WRK}/das/${NAME}/${NAME}_DASTool_bins 
NON_MOBILE_DAS=${WRK}/das/${NAME}/${NAME}_DASTool_non_mobile_bins

#run the python script to make this new table
python ~/agk/CDC2/scripts/combo_table_clusters.py -i $INFILE -o $OUTFILE -d $DAS -n $NON_MOBILE_DAS

