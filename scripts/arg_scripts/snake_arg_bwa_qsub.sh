#$ -S /bin/bash
#$ -N snk_bwa_arg
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/CDC2/logs/snake_bwa_arg_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/snake_bwa_arg_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2/scripts
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito2

FOLDER=CDC2
WRK="/workdir/users/agk85/${FOLDER}"
JOBS=44
RESTARTS=2
ERR="/workdir/users/agk85/${FOLDER}/logs/"
OUT=$ERR
#snakemake -s ${WRK}/scripts/snake_arg_index --jobs $JOBS --restart-times $RESTARTS --cluster "qsub -q long.q@cbsubrito2 -S /bin/bash -e $ERR -o $OUT -wd ${WRK}/arg_v_org/metagenomes/ -N {params.n} -l h_vem={resources.mem_mb}G"

#just run this to make sure you've indexed properly
#snakemake -s ${WRK}/scripts/arg_scripts/snake_arg_index

#should be quick!
snakemake -s ${WRK}/scripts/arg_scripts/snake_arg_bwa --jobs $JOBS --restart-times $RESTARTS --cluster "qsub -q long.q@cbsubrito2 -S /bin/bash -e $ERR -o $OUT -wd ${WRK}/arg_v_org/metagenomes/ -N {params.n} -l h_vmem={resources.mem_mb}G"
