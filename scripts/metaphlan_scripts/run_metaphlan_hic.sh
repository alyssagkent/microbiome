#$ -S /bin/bash
#$ -N mtphln_hic
#$ -V
#$ -e /workdir/users/agk85/CDC2/metaphlan/err/
#$ -o /workdir/users/agk85/CDC2/metaphlan/err/
#$ -t 5
#$ -wd /workdir/users/agk85/CDC2/metaphlan/
#$ -l h_vmem=40G
#$ -q long.q@cbsubrito2
#$ -pe parenv 1

## This script will run MetaPhlAn-2.0 to identify the organismal composition of 
##the CDC metagenomic data
#will need to redo 1 an 5 because B314-1 and B316-1

## Set directories
WRK=/workdir/users/agk85/CDC2
DATA=/workdir/data/CDC/hic/merged
OUT=$WRK/metaphlan/hic
REF=/workdir/users/fnn3/references

## Create design file of file names
LIST=$WRK/HicDesign.txt
NAME=$(sed -n "${SGE_TASK_ID}p" $LIST)

FASTQ1=$DATA/${NAME}hic.1.fastq
FASTQ2=$DATA/${NAME}hic.1.solo.fastq
FASTQ3=$DATA/${NAME}hic.2.fastq
FASTQ4=$DATA/${NAME}hic.2.solo.fastq

## Run metaphlan
export PATH=/programs/MetaPhlAn-2.0:/programs/MetaPhlAn-2.0/utils:$PATH

metaphlan2.py --input_type fastq --bowtie2db $REF/db_v20/mpa_v20_m200 --bowtie2out $OUT/${NAME}hic.bowtie2.bz2 --nproc 4 $FASTQ1,$FASTQ3 $OUT/${NAME}hic_profile.txt
