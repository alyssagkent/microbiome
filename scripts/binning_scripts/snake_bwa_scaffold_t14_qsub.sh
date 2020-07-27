#$ -S /bin/bash
#$ -N snake_bwa_t14
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake-bwa_t14_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake-bwa_t14_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/concoct
#$ -l h_vmem=1G
#$ -q long.q@cbsubrito2

#--restart-times 2
FOLDER=CDC2
SCRIPT=/workdir/users/agk85/${FOLDER}/scripts/binning_scripts/snake_bwa_scaffold_t14
RESTARTS=1
JOBS=5
ERR=/workdir/users/agk85/${FOLDER}/logs

snakemake -s $SCRIPT --restart-times $RESTARTS --jobs $JOBS --cluster "qsub -q long.q@cbsubrito2 -S /bin/bash -e $ERR -o $ERR -N {params.n} -l h_vmem={resources.mem_mb}G"
