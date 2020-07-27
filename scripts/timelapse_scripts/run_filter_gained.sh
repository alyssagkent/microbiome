
GENE=arg

python filter_gained.py /workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained.txt \
/workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained_filtered.txt \
/workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained_filtered_bestconnections.txt \
/workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained_filteredindex.txt

GENE=mge
python filter_gained.py /workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained.txt \
/workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained_filtered.txt \
/workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained_filtered_bestconnections.txt \
/workdir/users/agk85/CDC2/bins/timelapse/timelapse_${GENE}_org_2_alltaxalevels_gained_filteredindex.txt
