#$ -S /bin/bash
#$ -N plasflow_overlap
#$ -e /workdir/users/agk85/CDC2/logs/plasflow_overlap_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/plasflow_overlap_$JOB_ID.out
#$ -t 1-32
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=1G
#$ -q long.q@cbsubrito2

# Goal is to bin contigs into single species from idba generated scaffolds

WRK=/workdir/users/agk85/CDC2

DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

#run the overlap script on each file
python /workdir/users/agk85/CDC2/scripts/get_plasflow_overlap.py ${WRK}/combo_tables/metagenomes/${NAME}_master_scf_table.txt ${WRK}/combo_tables/metagenomes/${NAME}_plasflow_overlap.txt ${WRK}/combo_tables/metagenomes/${NAME}_plasmid_lengths.txt $NAME


