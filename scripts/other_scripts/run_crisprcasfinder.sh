#bash script to run crisprcasfinder

JOB_ID='manual'
for SGE_TASK_ID in {1..32};
do
WRK=/workdir/users/agk85/CDC2
WRK2=/workdir/agk85/crisprcasfinder
if [ ! -d $WRK2 ]; then mkdir -p $WRK2; fi

DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

cp $DESIGN_FILE $WRK2


OUT=$WRK2/${NAME};
if [ ! -d $OUT ]; then mkdir -p $OUT; fi

SCF=$WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta

cd $WRK2//${NAME}
cp $SCF $OUT

done



#This was to import the docker image?
#docker1 import /programs/CRISPRCasFinder/CRISPRCasFinder_20181004.tar
#docker1 run -d -it biohpc_agk85/crisprcasfinder_20181004 /bin/bash
#docker1 ps -a
#docker1 exec -it 51b8c42c907d /bin/bash
#source ~/.profile

#so once you are in the container your workdir is different
#WRK=/workdir/crisprcasfinder
#DESIGN_FILE=$WRK/MetaDesign.txt
#        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
#        NAME=`basename "$DESIGN"`
#SCF=$WRK/${NAME}/${NAME}_scaffold.fasta
#OUT=$WRK/${NAME}
#echo ${NAME};
#cd $OUT
#perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -out output -meta 1 -so  /usr/local/CRISPRCasFinder/sel392v2.so -i $SCF

#perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -o $OUT/temp -meta 1 -so  /usr/local/CRISPRCasFinder/sel392v2.so -i $SCF
#exit

#docker1 clean all


#######this is what idid with metabat for comparison
#docker1 run  metabat/metabat:latest metabat2 -i ${NAME}_scaffold.fasta -o ${NAME}_metabat -m 1500 -s 10000 --saveCls --unbinned --seed 23 ${NAME}.sorted.bam;
#docker1 clean


#something to move the unbinned to fasta so you can run DASTool without including the unbinned
#mv $WRK2/${NAME}_metabat* $WRK2/metabat/${NAME}
#mv $WRK2/metabat/${NAME}/${NAME}_metabat.unbinned.fa $WRK2/metabat/${NAME}/${NAME}_metabat.unbinned.fasta


#echo "docker1 run  metabat/metabat:latest metabat2 -i ${NAME}/${NAME}_scaffold.fasta -o ${NAME}/${NAME}_metabat -m 1500 -s 10000 --saveCls --unbinned --seed 23 ${NAME}/${NAME}.sorted.bam;" >> $LOG 2>&1
#echo end metabat ${NAME} `date` >> $LOG 2>&1

#done

