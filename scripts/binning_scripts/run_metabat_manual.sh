
#this script runs metabat on contigs
#this file requires mapped contigs for abundances
# it loops over the CDC and healthy patients
#then copies the scaffold
#then moves the sorted bam files mapped to the docker folder (/workdir/agk85/)
#then runs metabat
#and prints to a logfile in /workdir/users/agk85/CDC2/logs

JOB_ID='manual'
for SGE_TASK_ID in {7..16};
do
WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/MetaDesign_newmgms.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`


WRK2=/workdir/agk85
LOG=$WRK/logs/metabat_${NAME}_${JOB_ID}.out
echo ${NAME};
mkdir $WRK2/${NAME};
cp $WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta ${WRK2}/${NAME}/
cp $WRK/mapping/${NAME}/${NAME}.sorted.bam ${WRK2}/${NAME}/;
cp $WRK/mapping/${NAME}/${NAME}.sorted.bam.bai ${WRK2}/${NAME}/;

cd ${WRK2}

echo start metabat ${NAME} `date` >> $LOG 2>&1
docker1 run  metabat/metabat:latest metabat2 -i ${NAME}/${NAME}_scaffold.fasta -o ${NAME}/${NAME}_metabat -m 1500 -s 10000 --saveCls --unbinned --seed 23 ${NAME}/${NAME}.sorted.bam;
docker1 clean

#something to move the unbinned to fasta so you can run DASTool without including the unbinned
mv ${WRK2}/${NAME}/${NAME}_metabat.unbinned.fa ${WRK2}/${NAME}/${NAME}_metabat.unbinned.fasta
rm ${WRK2}/${NAME}/${NAME}.sorted.bam
rm ${WRK2}/${NAME}/${NAME}.sorted.bam.bai
rm ${WRK2}/${NAME}/${NAME}_scaffold.fasta
mv $WRK2/${NAME} $WRK2/metabat/${NAME}
echo end metabat ${NAME} `date` >> $LOG 2>&1
done

