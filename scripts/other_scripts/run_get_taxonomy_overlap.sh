#$ -S /bin/bash
#$ -N taxonomy_overlap
#$ -e /workdir/users/agk85/CDC2/logs/taxonomy_overlap_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/taxonomy_overlap_$JOB_ID.out
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
python /workdir/users/agk85/CDC2/scripts/get_taxonomy_overlap.py ${WRK}/combo_tables/metagenomes/${NAME}_master_scf_table.txt ${WRK}/combo_tables/metagenomes/${NAME}_taxonomy_overlap.txt $NAME


