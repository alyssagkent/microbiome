#$ -S /bin/bash
#$ -N rename
#$ -e /workdir/users/agk85/CDC2/logs/concoct_timepoints_rename_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/concoct_timepoints_rename_$JOB_ID.out
#$ -t 26-32
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 1
#$ -l h_vmem=1G
#$ -q short.q@cbsubrito2

#renames concoct_timepoints so that they have the sample name in them
WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

FOLDER=$WRK/concoct_timepoints/${NAME}/concoct_output/fasta_bins
cd $FOLDER
for file in *.fa;
do 
echo $file
echo ${NAME}_concoct_${file}
mv $file ${NAME}_concoct_${file}
done
