WRK=/workdir/users/agk85/CDC2
head -1 $WRK/bins/B314-1_bin_cleanliness.txt > $WRK/bins/all_bin_cleanliness_all.txt
tail -n +2 -q $WRK/bins/*_bin_cleanliness.txt >> $WRK/bins/all_bin_cleanliness_all.txt


Rscript $WRK/scripts/cleanliness_scripts/bin_cleanliness.R
