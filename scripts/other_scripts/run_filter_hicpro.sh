for SGE_TASK_ID in {1..44};
do
WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/HicDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`
echo ${NAME} `date` begin
INBAM=$WRK/hicpro/output/${NAME}_output/bowtie_results/bwt2/${NAME}/${NAME}_${NAME}_scaffold.fasta.bwt2pairs.bam
INSAM=$WRK/hicpro/output/${NAME}_output/bowtie_results/bwt2/${NAME}/${NAME}_${NAME}_scaffold.fasta.bwt2pairs.sam
OUTSAM_99=$WRK/hicpro/output/${NAME}_output/bowtie_results/bwt2/${NAME}/${NAME}_${NAME}_scaffold.fasta.bwt2pairs_99.sam
OUTSAM_98=$WRK/hicpro/output/${NAME}_output/bowtie_results/bwt2/${NAME}/${NAME}_${NAME}_scaffold.fasta.bwt2pairs_98.sam
#samtools view -h $INBAM > $INSAM
#perl /workdir/scripts/alignments/filterSamByIdentity.pl $INSAM 99 0 0 1 1 > $OUTSAM_99
python /workdir/users/agk85/CDC2/hicpro/scripts/which_trans_pid_trans.py $NAME 99

#perl /workdir/scripts/alignments/filterSamByIdentity.pl $INSAM 98 0 0 1 1 > $OUTSAM_98
python /workdir/users/agk85/CDC2/hicpro/scripts/which_trans_pid_trans.py $NAME 98

echo ${NAME} `date` done
done
