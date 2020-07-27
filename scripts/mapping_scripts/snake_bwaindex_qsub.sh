#$ -S /bin/bash
#$ -N snake_index
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake-bwaindex_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake-bwaindex_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/logs
#$ -l h_vmem=1G
#$ -q short.q@cbsubrito2

FOLDER=CDC2
SCRIPT=/workdir/users/agk85/${FOLDER}/scripts/snake_bwaindex
JOBS=8
RESTARTS=1
ERR=/workdir/users/agk85/${FOLDER}/logs


#--restart_times 2
snakemake -s $SCRIPT --jobs $JOBS --restart-times $RESTARTS --cluster "qsub -q short.q@cbsubrito2 -S /bin/bash -e $ERR -o $ERR -N {params.n} -l h_vmem={resources.mem_mb}G"
