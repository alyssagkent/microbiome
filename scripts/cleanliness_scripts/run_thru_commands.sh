
#this is running through org_v_org stuff to modify something in the annotation of things in the amphora combo_table section

#bash run_master_scf_table.sh
#bash run_master_scf_table_binary.sh
python /workdir/users/agk85/CDC2/scripts/cleanliness_scripts/hic_org_v_org.py /workdir/users/agk85/CDC2
Rscript /workdir/users/agk85/CDC2/scripts/cleanliness_scripts/org_v_org_CDC_cleanliness.R

#then go to the folder on computer /Users/agk/Box\ Sync/Labwork-CDC/hic/org_v_org
#scp agk85@cbsubrito2.tc.cornell.edu:~/agk/CDC2/figures/plotting/CDC_org_v_org* .


