WRK=/workdir/users/agk85/CDC2

CONTACT_THRESH=1
CONTACT_THRESHOLD=5
PATIENT=B314
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}


PATIENT=B316
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-6_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-7_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B320
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B331
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B335
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B357
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-6_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B370
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}


PATIENT=US3
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-8_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-10_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-12_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-14_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-16_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}


PATIENT=US8
python $WRK/scripts/mge_histogram_prep.py -i ${WRK}/bins/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD}


#################

PATIENT=all
python $WRK/scripts/mge_histogram_prep.py \
-i ${WRK}/bins/B314-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B314-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B314-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B314-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B316-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B316-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B316-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B316-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B316-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B316-6_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B316-7_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B320-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B320-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B320-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B320-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B331-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B331-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B331-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B331-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B335-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B335-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B335-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B357-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B357-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B357-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B357-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B357-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B357-6_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B370-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B370-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B370-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B370-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/B370-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US3-8_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US3-10_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US3-12_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US3-14_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US3-16_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US8-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US8-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US8-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US8-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins/US8-5_das_${CONTACT_THRESH}_contigs_v_bins.txt \
-a ${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl -p ${PATIENT} -c ${CONTACT_THRESHOLD} \
-o ${WRK}/bins/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt

#Rscript $WRK/scripts/histogram_mge_plotting.R $WRK/bins/${PATIENT}_das_${CONTACT_THRESH}_mgetaxa.txt $WRK/bins/${PATIENT}_das_${CONTACT_THRESH}_mgetaxa.txt 


cd  $WRK/bins

#head -1 Numbers_${CONTACT_THRESH}_B314.txt > ARG_histogram_numbers_${CONTACT_THRESH}.txt
#tail -n +2 -q Numbers*.txt >> ARG_histogram_numbers_${CONTACT_THRESH}.txt

#head -1 Taxa_Count_Numbers_${CONTACT_THRESH}_B314.txt > ARG_histogram_taxa_count_numbers_${CONTACT_THRESH}.txt
#tail -n +2 -q Taxa_Count_Numbers*.txt >> ARG_histogram_taxa_count_numbers_${CONTACT_THRESH}.txt

#head -1 Metrics_${CONTACT_THRESH}_B314.txt > ARG_histogram_metrics_${CONTACT_THRESH}.txt
#tail -n +2 -q Metrics*.txt >> ARG_histogram_metrics_${CONTACT_THRESH}.txt


head -1 ${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt > histogram_figures/Together_das_${CONTACT_THRESHOLD}_mgetaxa_together.txt
tail -n +2 -q *_das_${CONTACT_THRESHOLD}_mgetaxa.txt >> histogram_figures/Together_das_${CONTACT_THRESHOLD}_mgetaxa_together.txt



