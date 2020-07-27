#$ -S /bin/bash
#$ -N arg_cliques
#$ -e /workdir/users/agk85/CDC/arg_v_org/log/arg_cliques_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC/arg_v_org/log/arg_cliques_$JOB_ID.out
#$ -V
#$ -t 1-3
#$ -wd /workdir/users/agk85/CDC
#$ -pe parenv 1
#$ -l h_vmem=30G
#$ -q long.q@cbsubrito2

# Goal is to run the mge_org_hic_refined_patient_hicclusters_long.py file with the different paramters to get mge counts 
#run python script
echo $SGE_TASK_ID

depth=$(echo "$(($SGE_TASK_ID-1))")
SCRIPT=/workdir/users/agk85/CDC/arg_v_org/metagenomes3/scripts/arg_org_hic_refined_patient_hicclusters_cliques_v4.py
python $SCRIPT -d $depth

