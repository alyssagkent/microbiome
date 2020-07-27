#$ -S /bin/bash
#$ -N hicpro
#$ -V
#$ -t 33-44
#$ -e /workdir/users/agk85/CDC2/logs/hicpro_$JOB_ID.e
#$ -o /workdir/users/agk85/CDC2/logs/hicpro_$JOB_ID.e
#$ -wd /workdir/users/agk85/CDC2
#$ -l h_vmem=30G
#$ -q short.q@cbsubrito2

for SGE_TASK_ID in {33..44};
do

#This script runs hicpro setup.
WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/HicDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

SCF=$WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta
SCFNAME=${NAME}_scaffold.fasta
cd $WRK/hicpro/references

export PYTHONPATH=/programs/HiC-Pro_2.7.9/lib64/python2.7/site-packages:$PYTHONPATH
export PATH=/programs/HiC-Pro_2.7.9/bin:$PATH


#I did this manually
bowtie2-build $SCF $SCFNAME
python ~/agk/CDC/scripts/reference_genome_sizes.py $SCF ${NAME}_scaffold.sizes
python /programs/HiC-Pro_2.7.9/bin/utils/digest_genome.py -r dpnii -o ${NAME}_sau3ai.bed $SCF
OUTPUT=$WRK/hicpro/rawdata/${NAME}_rawdata/${NAME}
if [ ! -d $OUTPUT ]; then mkdir -p $OUTPUT; fi
#R1=/workdir/data/CDC/hic/merged/${NAME}hic.1.fastq
#R2=/workdir/data/CDC/hic/merged/${NAME}hic.2.fastq
#cp $R1 ${OUTPUT}/${NAME}_R1.fastq
#cp $R2 ${OUTPUT}/${NAME}_R2.fastq
#/programs/HiC-Pro_2.7.9/bin/utils/split_reads.py --results_folder $WRK/hicpro/${NAME}_rawdata/${NAME} --nreads 500000 $OUTPUT/${NAME}_R1.fastq
#/programs/HiC-Pro_2.7.9/bin/utils/split_reads.py --results_folder $WRK/hicpro/${NAME}_rawdata/${NAME} --nreads 500000 $OUTPUT/${NAME}_R2.fastq

#cp /workdir/users/agk85/CDC2/hicpro/configs/sau3ai_config-hicpro.txt /workdir/users/agk85/CDC2/hicpro/configs/${NAME}_config-hicpro.txt
#sed -i "s/B314-1/$NAME/g" "/workdir/users/agk85/CDC2/hicpro/configs/${NAME}_config-hicpro.txt"
done
