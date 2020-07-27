#$ -S /bin/bash
#$ -N snake_r10l1
#$ -V
#$ -t 1
#$ -e /workdir/users/agk85/cdc_QC/round10/lane1/err/snakemake_$JOB_ID.err
#$ -o /workdir/users/agk85/cdc_QC/round10/lane1/log/snakemake_$JOB_ID.out
#$ -wd /workdir/users/agk85/cdc_QC/round10/lane1
#$ -l h_vmem=1G
#$ -q long.q@cbsubrito2

snakemake -s /workdir/users/agk85/cdc_QC/round10/lane1/snake_QC --jobs 16 --restart-times 5 --cluster "qsub -q long.q@cbsubrito2 -S /bin/bash -e /workdir/users/agk85/cdc_QC/round10/lane1/err -o /workdir/users/agk85/cdc_QC/round10/lane1/log -N {params.n} -l h_vmem={resources.mem_mb}G"
