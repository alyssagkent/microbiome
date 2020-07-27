
cd ~/agk/CDC2/maxbin
for fold in B*; 
do 
awk -v OFS="\t" '$1=$1' ${fold}/checkm_lineage/${fold}.qa | grep -v '^Bin' | cut -f1,2,7,8,9,10 > ${fold}/${fold}.stats;
done

for fold in U*;
do
awk -v OFS="\t" '$1=$1' ${fold}/checkm_lineage/${fold}.qa | grep -v '^Bin' | cut -f1,2,7,8,9,10 > ${fold}/${fold}.stats;
done
