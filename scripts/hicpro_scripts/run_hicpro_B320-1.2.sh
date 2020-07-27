#$ -S /bin/bash
#$ -N hicpro
#$ -V
#$ -e /workdir/users/agk85/CDC2/logs/hicpro_pieces_$JOB_ID.out
#$ -o /workdir/users/agk85/CDC2/logs/hicpro_pieces_$JOB_ID.out
#$ -wd /workdir/users/agk85/CDC2
#$ -l h_vmem=200G
#$ -q long.q@cbsubrito2

WRK=/workdir/users/agk85/CDC2
NAME=B320-1.2

export PYTHONPATH=/programs/HiC-Pro_2.7.9/lib64/python2.7/site-packages:$PYTHONPATH
export PATH=/programs/HiC-Pro_2.7.9/bin:$PATH
echo started hic-pro ${NAME} `date`
INPUT=$WRK/hicpro/${NAME}_rawdata
OUTPUT=$WRK/hicpro/output/${NAME}_output
CONFIG=$WRK/hicpro/configs/${NAME}_config-hicpro.txt


#/programs/HiC-Pro_2.7.9/bin/utils/split_reads.py --results_folder $INPUT/${NAME} --nreads 1000000 $INPUT/${NAME}/${NAME}_R1.fastq
#/programs/HiC-Pro_2.7.9/bin/utils/split_reads.py --results_folder $INPUT/${NAME} --nreads 1000000 $INPUT/${NAME}/${NAME}_R2.fastq

echo started hic-pro `date`
#Run hic-pro
#HiC-Pro -i $INPUT -o $OUTPUT -c $CONFIG

#stepwise
HiC-Pro -i $INPUT -o $OUTPUT -c $CONFIG -s mapping -s quality_checks
echo finished hic-pro mapping $NAME `date`
HiC-Pro -i $OUTPUT/bowtie_results/bwt2 -o $OUTPUT -c $CONFIG -s proc_hic -s quality_checks
echo finished hic-pro proc_hic $NAME `date`
HiC-Pro -i $OUTPUT/hic_results/data -o $OUTPUT -c $CONFIG -s merge_persample -s quality_checks
echo finished hic-pro merge `date`
