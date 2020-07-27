#$ -S /bin/bash
#$ -N trim_pro_phage_is
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake_min_prodigal_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake_min_prodigal_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/logs/
#$ -l h_vmem=1G
#$ -q short.q@cbsubrito2

#--restart-times 2

FOLDER=CDC2
SCRIPT=/workdir/users/agk85/${FOLDER}/scripts/snakes/snake_min_prodigal_pf
JOBS=5
RESTARTS=2
LOG=/workdir/users/agk85/${FOLDER}/logs
snakemake -s $SCRIPT --jobs $JOBS --restart-times $RESTARTS --cluster "qsub -q short.q@cbsubrito2 -S /bin/bash -e $LOG -o $LOG -N {params.n} -l h_vmem={resources.mem_mb}G"

