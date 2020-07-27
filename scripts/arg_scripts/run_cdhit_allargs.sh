#$ -S /bin/bash
#$ -N cdhit
#$ -e /workdir/users/agk85/CDC2/logs/cdhit_arg_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/cdhit_arg_$JOB_ID.out
#$ -V
#$ -pe parenv 1
#$ -l h_vmem=30G
#$ -t 1
#$ -q short.q@cbsubrito2

# Program will cluster all ARGs from a patient at 99% identity using cdhit 
echo start cdhit `date`
FOLDER=CDC2
NUM=99
PERCENTAGE=0.${NUM}
#change this if it is 100 because 0.1 is wrong hah
#PERCENTAGE=.9
echo $NUM
#run cd-hit using output from prodigal that has been concatenated to include all metagenomes ORFs by patient in respective files
/programs/cd-hit-v4.6.1-2012-08-27/cd-hit-est -i /workdir/users/agk85/${FOLDER}/args/args.fasta -o /workdir/users/agk85/${FOLDER}/args/args_${NUM}_nr.fna -c $PERCENTAGE -M 30000 -T 8 -n 8 -d 0 -s .9

#make it into a nice table
python $WRK/scripts/arg_scripts/cdhit_output_to_table.py /workdir/users/agk85/CDC2/args/args_99_nr.fna.clstr
echo end cd-hit `date`
