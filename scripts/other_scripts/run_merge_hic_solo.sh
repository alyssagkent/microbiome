#$ -S /bin/bash
#$ -N merge_hic_solo
#$ -V
#$ -t 4
#$ -e /workdir/data/CDC/hic/merged/err
#$ -o /workdir/data/CDC/hic/merged/log
#$ -wd /workdir/users/agk85/CDC
#$ -l h_vmem=20G
#$ -q long.q@cbsubrito2


#This script merges gzipped files

WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/HicDesign_new.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

cat /home/britolab/data/CDC/trimmedReads/round*/hic/lane*/solo/${NAME}hic.1.solo.fastq.gz > /workdir/data/CDC/hic/merged/${NAME}hic.1.solo.fastq.gz
cat /home/britolab/data/CDC/trimmedReads/round*/hic/lane*/solo/${NAME}hic.2.solo.fastq.gz > /workdir/data/CDC/hic/merged/${NAME}hic.2.solo.fastq.gz
gunzip /workdir/data/CDC/hic/merged/${NAME}hic.1.solo.fastq.gz
gunzip /workdir/data/CDC/hic/merged/${NAME}hic.2.solo.fastq.gz

