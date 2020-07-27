#$ -S /bin/bash
#$ -N cdhit
#$ -e /workdir/users/agk85/CDC/arg_v_org/log/cdhit_arg_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC/arg_v_org/log/cdhit_arg_$JOB_ID.out
#$ -V
#$ -pe parenv 8
#$ -l h_vmem=30G
#$ -t 1
#$ -q short.q@cbsubrito2

# Program will cluster all ARGs from a patient at 99% identity using cdhit 
echo start cdhit `date`
NUM=100
PERCENTAGE=0.${NUM}
PERCENTAGE=1.0
echo $NUM
#run cd-hit using output from prodigal that has been concatenated to include all metagenomes ORFs by patient in respective files
/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i /workdir/users/agk85/CDC/arg_v_org/metagenomes3/arg_prot.fasta -o /workdir/users/agk85/CDC/arg_v_org/metagenomes3/args_${NUM}_nr.fna -c $PERCENTAGE -M 30000 -T 8 -n 8 -d 0 -s .9

echo end cd-hit `date`
