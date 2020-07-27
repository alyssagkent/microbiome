#$ -S /bin/bash
#$ -N metagenome_stats
#$ -V
#$ -t 1-43
#$ -e /workdir/users/agk85/CDC2/logs/quasts_${JOB_ID}.err
#$ -wd /workdir/users/agk85/CDC
#$ -l h_vmem=5G
#$ -pe parenv 8
#$ -q short.q@cbsubrito2



#Run QUAST to get assembly statistics

WRK=/workdir/users/agk85


LIST=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $LIST)
        NAME=`basename "$DESIGN"`

ASSEMBLY=$WRK/idba_rerun/metagenomes/${NAME}/${NAME}_scaffolds.fasta
OUT=$WRK/idba_rerun/metagenomes/${NAME}/${NAME}_assembly_stats
        if [ ! -d $OUT ]; then mkdir -p $OUT; fi

python /programs/quast-4.0/quast.py -o $OUT $ASSEMBLY
