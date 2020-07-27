# Identify mobile taxonomic associations

echo start blast `date`

WRK=/workdir/users/agk85/CDC2

OUT=$WRK/mobile/metagenomes/${NAME}
if [ ! -d $OUT ]; then mkdir -p $OUT; fi

DB=/workdir/blastdb/nt

#made this with commandline input
#makeblastdb -in /workdir/refdbs/plasmids/plasmids.fna -dbtype 'nucl' -parse_seqids -out plasmids_db
QUERY=$WRK/mobile/metagenomes/mge_99_nr.fna

#run blastp using output from prodigal that has been concatenated to include all metagenomes and all phage in respective files
/programs/ncbi-blast-2.3.0+/bin/blastn -query $QUERY -db $DB -out $OUT/mge_99_nr.out -num_threads 22 -outfmt "6 qseqid sseqid stitle sblastnames salltitles pident length mismatch staxids" -evalue 1e-10

