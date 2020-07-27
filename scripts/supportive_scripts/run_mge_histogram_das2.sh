WRK=/workdir/users/agk85/CDC2
CONTACT_THRESH=1
CONTACT_THRESHOLD=2
SCRIPT=$WRK/scripts/supportive_scripts/mge_histogram_prep.py
CLUSTER=${WRK}/mobile/metagenomes/mge_99_nr.fna.clstr.tbl


PATIENT=B314
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}


PATIENT=B316
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-6_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-7_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B320
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B331
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B335
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B357
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-6_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=B370
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}


PATIENT=US3
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-8_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-10_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-12_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-14_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-16_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}

PATIENT=US8
python $SCRIPT -i ${WRK}/bins_hicsupport/${PATIENT}-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/${PATIENT}-5_das_${CONTACT_THRESH}_contigs_v_bins.txt -o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt -a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD}


#################

PATIENT=all
python $SCRIPT \
-i ${WRK}/bins_hicsupport/B314-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B314-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B314-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B314-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B316-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B316-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B316-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B316-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B316-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B316-6_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B316-7_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B320-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B320-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B320-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B320-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B331-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B331-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B331-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B331-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B335-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B335-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B335-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B357-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B357-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B357-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B357-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B357-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B357-6_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B370-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B370-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B370-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B370-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/B370-5_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US3-8_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US3-10_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US3-12_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US3-14_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US3-16_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US8-1_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US8-2_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US8-3_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US8-4_das_${CONTACT_THRESH}_contigs_v_bins.txt,${WRK}/bins_hicsupport/US8-5_das_${CONTACT_THRESH}_contigs_v_bins.txt \
-a $CLUSTER -p ${PATIENT} -c ${CONTACT_THRESHOLD} \
-o ${WRK}/bins_hicsupport/${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt

cd  $WRK/bins_hicsupport

head -1 ${PATIENT}_das_${CONTACT_THRESHOLD}_mgetaxa.txt > Together_das_${CONTACT_THRESHOLD}_mgetaxa_together.txt
tail -n +2 -q *_das_${CONTACT_THRESHOLD}_mgetaxa.txt >> Together_das_${CONTACT_THRESHOLD}_mgetaxa_together.txt

head -1 B314-1_das_1_contigs_v_bins.txt > all_das_1_contigs_v_bins_all.txt
tail -n +2 -q *_das_1_contigs_v_bins.txt >> all_das_1_contigs_v_bins_all.txt


