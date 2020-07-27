


WRK=/workdir/users/agk85/CDC2
bash $WRK/scripts/supportive_scripts/run_arg_histogram_das2.sh
bash $WRK/scripts/supportive_scripts/run_arg_histogram_das5.sh
bash $WRK/scripts/supportive_scripts/run_mge_histogram_das2.sh
bash $WRK/scripts/supportive_scripts/run_mge_histogram_das5.sh

python $WRK/scripts/supportive_scripts/linelists.py

cd /workdir/users/agk85/CDC2/bins_hicsupport
rm Together_das_*linelist.txt
cat *_das_2_argtaxa_linelist.txt > Together_das_2_argtaxa_linelist.txt
cat *_das_5_argtaxa_linelist.txt > Together_das_5_argtaxa_linelist.txt
cat *_das_2_mgetaxa_linelist.txt > Together_das_2_mgetaxa_linelist.txt
cat *_das_5_mgetaxa_linelist.txt > Together_das_5_mgetaxa_linelist.txt


Rscript ~/agk/CDC2/scripts/supportive_scripts/plotting_supportive.R
