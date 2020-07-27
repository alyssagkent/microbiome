#$ -S /bin/bash
#$ -N plasflow
#$ -e /workdir/users/agk85/CDC2/logs/plasflow_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/plasflow_$JOB_ID.out
#$ -t 1-37
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito2

# Goal is to bin contigs into single species from idba generated scaffolds

WRK=/workdir/users/agk85/CDC2

DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`


/programs/mauve_2.4.0/linux-x64/progressiveMauve --output=/workdir/users/agk85/CDC/resfams/metagenomes/mauve/${clusternum} $file
