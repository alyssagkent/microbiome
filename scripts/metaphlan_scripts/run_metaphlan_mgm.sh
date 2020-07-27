#$ -S /bin/bash
#$ -N CDC2_phln
#$ -V
#$ -e /workdir/users/agk85/CDC2/metaphlan/err/
#$ -o /workdir/users/agk85/CDC2/metaphlan/err/
#$ -t 33-44
#$ -wd /workdir/users/agk85/CDC2/metaphlan/
#$ -l h_vmem=35G
#$ -q short.q@cbsubrito2
#$ -pe parenv 4

## This script will run MetaPhlAn-2.0 to identify the organismal composition of 
##the CDC metagenomic data

#tasks to do B320-2 = 13, B331-2=18

## Set directories
WRK=/workdir/users/agk85/CDC2
DATA=/workdir/data/CDC/metagenomes/merged
OUT=$WRK/metaphlan/mgm
REF=/workdir/users/fnn3/references
###set to HicDesign.txt to just do those ones quickly

## Create design file of file names
LIST=$WRK/MetaDesign.txt
NAME=$(sed -n "${SGE_TASK_ID}p" $LIST)

FASTQ1=$DATA/unzip/${NAME}.1.fastq
FASTQ2=$DATA/solo/${NAME}.1.solo.fastq
FASTQ3=$DATA/unzip/${NAME}.2.fastq
FASTQ4=$DATA/solo/${NAME}.2.solo.fastq

## Run metaphlan
export PATH=/programs/MetaPhlAn-2.0:/programs/MetaPhlAn-2.0/utils:$PATH
metaphlan2.py --input_type fastq --bowtie2db $REF/db_v20/mpa_v20_m200 --bowtie2out $OUT/${NAME}.bowtie2.bz2 --nproc 4 $FASTQ1,$FASTQ2,$FASTQ3,$FASTQ4 $OUT/${NAME}_profile.txt

