#$ -S /bin/bash
#$ -N snake_metaspades
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake_metaspades_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake_metaspades_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/scripts
#$ -l h_vmem=1G
#$ -q long.q@cbsubrito2

FOLDER=CDC2
SCRIPT=/workdir/users/agk85/${FOLDER}/scripts/snake_metaspades
JOBS=10
RESTARTS=2
LOG=/workdir/users/agk85/${FOLDER}/logs
snakemake -s $SCRIPT --jobs $JOBS --restart-times $RESTARTS --cluster "qsub -q long.q@cbsubrito2 -pe parenv 16 -S /bin/bash -e $LOG -o $LOG -N {params.n} -l h_vmem={resources.mem_mb}G"
