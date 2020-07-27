
#run some things before cd-hit in args
FOLDER=CDC2
WRK=/workdir/users/agk85/${FOLDER}

mkdir ${WRK}/args
#get card ids
cat ${WRK}/card/metagenomes/*/*.txt | cut -f1,9 | sort | uniq | grep -v '^ORF' > ${WRK}/args/card_prot_id_type.txt
cat ${WRK}/card/metagenomes/*/*.txt | cut -f1 | sort | uniq | grep -v '^ORF'  > ${WRK}/args/card_protein_ids.txt
#because CARD adds a stupid space after their ids
sed -i -e 's/ \t/\t/g' ${WRK}/args/card_prot_id_type.txt
#get resfams ids
cat ${WRK}/resfams/metagenomes/*/*_resfams.txt | cut -f1,3 | sort | uniq > ${WRK}/args/resfams_prot_id_type.txt
cat ${WRK}/resfams/metagenomes/*/*_resfams.txt | cut -f1 | sort | uniq > ${WRK}/args/resfams_protein_ids.txt

#combine them and get rid of the space after the CARD ids
#because CARD adds a stupid space after their ids
sed -i -e 's/ \t/\t/g' ${WRK}/args/card_prot_id_type.txt
cat ${WRK}/args/card_protein_ids.txt ${WRK}/args/resfams_protein_ids.txt | tr -d " \t\r" | sort | uniq > ${WRK}/args/arg_protein_ids.txt

#get all of the proteins in one reference file
cat ${WRK}/prodigal_excise/metagenomes/*/*_proteins.fna > ${WRK}/prodigal_excise/metagenomes/all_proteins.fna
python /workdir/users/agk85/CDC/scripts/general_getseqs.py ${WRK}/prodigal_excise/metagenomes/all_proteins.fna ${WRK}/args/arg_protein_ids.txt ${WRK}/args/args.fasta 0 all



#cluster them update this so it reflects the correct directory and identity threshold
#qsub ${WRK}/scripts/arg_scripts/run_cdhit_allargs.sh
#map them, update this

#ok going to try indexing within the bwa qsub script
#qsub ${WRK}/scripts/arg_scripts/snake_arg_bwa_qsub.sh

#get the arg names and make a 80% coverage table of the rpkms
#python ${WRK}/scripts/arg_scripts/arg_sample_rpkm.py ${WRK}/args/mapping/bwa_alignments_99_99 ${WRK}/args/mapping/bwa_alignments_99_99/arg_v_samp_99_99.txt ${WRK}/args/mapping/bwa_alignments_99_99/arg_v_samp_99_99_names.txt ${FOLDER}


