#$ -S /bin/bash
#$ -N Concoct_t5
#$ -e /workdir/users/agk85/CDC2/logs/concoct_t5_$JOB_ID.err
#$ -o /workdir/users/agk85/CDC2/logs/concoct_t5_$JOB_ID.out
#$ -t 1-4
#$ -V
#$ -wd /workdir/users/agk85/CDC2
#$ -pe parenv 8
#$ -l h_vmem=16G
#$ -q long.q@cbsubrito2

#This runs Concoct to bin contigs

#files
WRK=/workdir/users/agk85/CDC2
DESIGN_FILE=$WRK/concoct/MetaDesign_t5.txt
        DESIGN=$(sed -n "${SGE_TASK_ID}p" $DESIGN_FILE)
        NAME=`basename "$DESIGN"`

OUT=$WRK/concoct_timepoints/${NAME}
 if [ ! -d $OUT ]; then mkdir -p $OUT; fi
PATIENT=$(echo $NAME | cut -d'-' -f1);

export PATH=/programs/Anaconda2/bin:$PATH
export LD_LIBRARY_PATH=/programs/Anaconda2/lib:$LD_LIBRARY_PATH

#prep the scaffolds
/workdir/users/agk85/tools/CONCOCT/scripts/cut_up_fasta.py $WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta -c 10000 -o 0 --merge_last -b ${OUT}/${NAME}_contigs_10K.bed > ${OUT}/${NAME}_contigs_10K.fa

source activate test-samtools

#using * to get all of the folders correctamente!
/workdir/users/agk85/tools/CONCOCT/scripts/concoct_coverage_table.py ${OUT}/${NAME}_contigs_10K.bed $WRK/mapping/mapping_t5/${PATIENT}*/*.sorted.bam  > ${OUT}/${NAME}_coverage_table.tsv

conda deactivate

source activate concoct_env
concoct --composition_file ${OUT}/${NAME}_contigs_10K.fa --coverage_file ${OUT}/${NAME}_coverage_table.tsv -b ${OUT}/concoct_output/

conda deactivate 

/workdir/users/agk85/tools/CONCOCT/scripts/merge_cutup_clustering.py ${OUT}/concoct_output/clustering_gt1000.csv > ${OUT}/concoct_output/clustering_merged.csv

mkdir ${OUT}/concoct_output/fasta_bins

/workdir/users/agk85/tools/CONCOCT/scripts/extract_fasta_bins.py $WRK/prodigal_excise/metagenomes/${NAME}/${NAME}_scaffold.fasta ${OUT}/concoct_output/clustering_merged.csv --output_path ${OUT}/concoct_output/fasta_bins



