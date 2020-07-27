#$ -S /bin/bash
#$ -N arg_figures
#$ -e /workdir/users/agk85/CDC/arg_v_org/log/arg_figures_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC/arg_v_org/log/arg_figures_$JOB_ID.out
#$ -t 3
#$ -V
#$ -wd /workdir/users/agk85/CDC
#$ -l h_vmem=10G
#$ -q long.q@cbsubrito2

# Goal is to filter argannot output 100% identity and 100% coverage and convert to gff3 file

WRK=/workdir/users/agk85/CDC
DESIGN_FILE=$WRK/arg_v_org/metagenomes3/scripts/scriptfiles.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        SCRIPT=`basename "$DESIGN"`

echo $SCRIPT
Rscript ${WRK}/arg_v_org/metagenomes3/scripts/$SCRIPT

