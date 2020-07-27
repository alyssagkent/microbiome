#$ -S /bin/bash
#$ -N snake_id
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake_identify_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake_identify_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/prodigal_excise
#$ -l h_vmem=1G
#$ -q long.q@cbsubrito2

FOLDER=CDC2
SCRIPT=/workdir/users/agk85/${FOLDER}/scripts/snakes/snake_identify
JOBS=10
RESTARTS=1
LOG=/workdir/users/agk85/${FOLDER}/logs
snakemake -s $SCRIPT --jobs $JOBS --restart-times $RESTARTS --cluster "qsub -q long.q@cbsubrito2 -pe parenv {params.j} -S /bin/bash -e $LOG -o $LOG -N {params.n} -l h_vmem={resources.mem_mb}G"
