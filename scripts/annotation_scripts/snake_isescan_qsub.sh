#$ -S /bin/bash
#$ -N snake_isescan
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake_isescan_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake_isescan_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/iselements
#$ -l h_vmem=1G
#$ -q short.q@cbsubrito2

FOLDER=CDC2
SCRIPT=/workdir/users/agk85/${FOLDER}/scripts/snakes/snake_isescan
JOBS=11
RESTARTS=1
LOG=/workdir/users/agk85/${FOLDER}/logs
snakemake -s $SCRIPT --jobs $JOBS --restart-times $RESTARTS --cluster "qsub -q long.q@cbsubrito2 -pe parenv {params.j} -S /bin/bash -e $LOG -o $LOG -N {params.n} -l h_vmem={resources.mem_mb}G"
