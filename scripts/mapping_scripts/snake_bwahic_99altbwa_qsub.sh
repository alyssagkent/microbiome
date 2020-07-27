#$ -S /bin/bash
#$ -N snake_bwahic
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake-bwahic_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake-bwahic_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/scripts
#$ -l h_vmem=1G
#$ -q long.q@cbsubrito2

#--restart-times 2
FOLDER=CDC2
SCRIPT=/workdir/users/agk85/${FOLDER}/scripts/snake_bwahic_99altbwa
RESTARTS=1
JOBS=32
ERR=/workdir/users/agk85/${FOLDER}/logs

snakemake -s $SCRIPT --restart-times $RESTARTS --jobs $JOBS --cluster "qsub -q long.q@cbsubrito2 -S /bin/bash -e $ERR -o $ERR -N {params.n} -l h_vmem={resources.mem_mb}G"
