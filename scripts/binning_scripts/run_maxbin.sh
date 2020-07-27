#$ -S /bin/bash
#$ -N maxbin
#$ -e /workdir/users/agk85/CDC2/logs/maxbin_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/maxbin_$JOB_ID.out
#$ -t 23 
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 30
#$ -l h_vmem=50G
#$ -q short.q@cbsubrito2

# Goal is to bin contigs into single species from idba generated scaffolds

WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/MetaDesign.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

OUT=$WRK/maxbin/${NAME}
 if [ ! -d $OUT ]; then mkdir -p $OUT; fi

LOG=$WRK/logs/maxbin_$JOB_ID_${NAME}.out

SCF=$WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta
R1=/workdir/data/CDC/metagenomes/merged/unzip/${NAME}.1.fastq
R2=/workdir/data/CDC/metagenomes/merged/unzip/${NAME}.2.fastq
echo start maxbin ${NAME} `date` >> $LOG 2>&1
export PATH=/programs/bowtie2-2.3.0:$PATH
export PATH=/programs/idba-1.1.1/bin:$PATH
export PATH=/programs/hmmer/binaries:$PATH
#run maxbin
perl /programs/MaxBin-2.2.4/run_MaxBin.pl -contig $SCF -reads $R1 -reads2 $R2 -max_iteration 50 -thread 8 -out $OUT/${NAME} >> $LOG 2>&1

echo end maxbin ${NAME} `date` >> $LOG 2>&1

