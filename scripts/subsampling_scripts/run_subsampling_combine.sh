#!/bin/bash

SCRIPT=/workdir/users/agk85/CDC2/scripts/subsampling_scripts/subsampling_plotting.R
cd /workdir/users/agk85/CDC2/resampling

head -1 B314-1_subsampling_2_arg_conbin.txt > all_2_conbin.txt
tail -n +2 -q *_subsampling_2_arg_conbin.txt >> all_2_conbin.txt

head -1 B314-1_subsampling_2_arg_conbin.txt > all_5_conbin.txt
tail -n +2 -q *_subsampling_5_arg_conbin.txt >> all_5_conbin.txt

head -1 B314-1_subsampling_2_arg_connections.txt > all_2_arg_connections.txt
tail -n +2 -q *_subsampling_2_arg_connections.txt >> all_2_arg_connections.txt

head -1 B314-1_subsampling_5_arg_connections.txt > all_5_arg_connections.txt
tail -n +2 -q *_subsampling_5_arg_connections.txt >> all_5_arg_connections.txt

head -1 B314-1_subsampling_2_mge_connections.txt > all_2_mge_connections.txt
tail -n +2 -q *_subsampling_2_mge_connections.txt >> all_2_mge_connections.txt

head -1 B314-1_subsampling_5_mge_connections.txt > all_5_mge_connections.txt
tail -n +2 -q *_subsampling_5_mge_connections.txt >> all_5_mge_connections.txt


head -1 B314-1_subsampling_2_arg_conbin_vp.txt > all_2_conbin_vp.txt
tail -n +2 -q *_subsampling_2_arg_conbin_vp.txt >> all_2_conbin_vp.txt

head -1 B314-1_subsampling_2_arg_conbin_vp.txt > all_5_conbin_vp.txt
tail -n +2 -q *_subsampling_5_arg_conbin_vp.txt >> all_5_conbin_vp.txt

head -1 B314-1_subsampling_2_arg_connections_vp.txt > all_2_arg_connections_vp.txt
tail -n +2 -q *_subsampling_2_arg_connections_vp.txt >> all_2_arg_connections_vp.txt

head -1 B314-1_subsampling_5_arg_connections_vp.txt > all_5_arg_connections_vp.txt
tail -n +2 -q *_subsampling_5_arg_connections_vp.txt >> all_5_arg_connections_vp.txt

head -1 B314-1_subsampling_2_mge_connections_vp.txt > all_2_mge_connections_vp.txt
tail -n +2 -q *_subsampling_2_mge_connections_vp.txt >> all_2_mge_connections_vp.txt

head -1 B314-1_subsampling_5_mge_connections_vp.txt > all_5_mge_connections_vp.txt
tail -n +2 -q *_subsampling_5_mge_connections_vp.txt >> all_5_mge_connections_vp.txt


head -1 B314-1_subsampling_2_arg_conbin_trimmed.txt > all_2_conbin_trimmed.txt
tail -n +2 -q *_subsampling_2_arg_conbin_trimmed.txt >> all_2_conbin_trimmed.txt

head -1 B314-1_subsampling_2_arg_conbin_trimmed.txt > all_5_conbin_trimmed.txt
tail -n +2 -q *_subsampling_5_arg_conbin_trimmed.txt >> all_5_conbin_trimmed.txt

head -1 B314-1_subsampling_2_arg_connections_trimmed.txt > all_2_arg_connections_trimmed.txt
tail -n +2 -q *_subsampling_2_arg_connections_trimmed.txt >> all_2_arg_connections_trimmed.txt

head -1 B314-1_subsampling_5_arg_connections_trimmed.txt > all_5_arg_connections_trimmed.txt
tail -n +2 -q *_subsampling_5_arg_connections_trimmed.txt >> all_5_arg_connections_trimmed.txt

head -1 B314-1_subsampling_2_mge_connections_trimmed.txt > all_2_mge_connections_trimmed.txt
tail -n +2 -q *_subsampling_2_mge_connections_trimmed.txt >> all_2_mge_connections_trimmed.txt

head -1 B314-1_subsampling_5_mge_connections_trimmed.txt > all_5_mge_connections_trimmed.txt
tail -n +2 -q *_subsampling_5_mge_connections_trimmed.txt >> all_5_mge_connections_trimmed.txt




Rscript $SCRIPT all_2_arg_connections.txt all_2_arg_connections.pdf
Rscript $SCRIPT all_5_arg_connections.txt all_5_arg_connections.pdf
Rscript $SCRIPT all_2_mge_connections.txt all_2_mge_connections.pdf
Rscript $SCRIPT all_5_mge_connections.txt all_5_mge_connections.pdf
Rscript $SCRIPT all_2_conbin.txt all_2_conbin.pdf
Rscript $SCRIPT all_5_conbin.txt all_5_conbin.pdf


Rscript $SCRIPT all_2_arg_connections_vp.txt all_2_arg_connections_vp.pdf
Rscript $SCRIPT all_5_arg_connections_vp.txt all_5_arg_connections_vp.pdf
Rscript $SCRIPT all_2_mge_connections_vp.txt all_2_mge_connections_vp.pdf
Rscript $SCRIPT all_5_mge_connections_vp.txt all_5_mge_connections_vp.pdf
Rscript $SCRIPT all_2_conbin_vp.txt all_2_conbin_vp.pdf
Rscript $SCRIPT all_5_conbin_vp.txt all_5_conbin_vp.pdf

Rscript $SCRIPT all_2_arg_connections_trimmed.txt all_2_arg_connections_trimmed.pdf
Rscript $SCRIPT all_5_arg_connections_trimmed.txt all_5_arg_connections_trimmed.pdf
Rscript $SCRIPT all_2_mge_connections_trimmed.txt all_2_mge_connections_trimmed.pdf
Rscript $SCRIPT all_5_mge_connections_trimmed.txt all_5_mge_connections_trimmed.pdf
Rscript $SCRIPT all_2_conbin_trimmed.txt all_2_conbin_trimmed.pdf
Rscript $SCRIPT all_5_conbin_trimmed.txt all_5_conbin_trimmed.pdf
