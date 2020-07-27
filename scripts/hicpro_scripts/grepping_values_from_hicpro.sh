grep 'processed' *_output/bowtie_results/bwt2/*/*mpairstat
grep 'Unmapped_pairs' *_output/bowtie_results/bwt2/*/*mpairstat
grep 'Unique_paired_alignments' *_output/bowtie_results/bwt2/*/*mpairstat
grep 'Multiple_pairs_alignments' *_output/bowtie_results/bwt2/*/*mpairstat
grep 'Pairs_with_singleton' *_output/bowtie_results/bwt2/*/*mpairstat

grep 'valid_interaction'$'\t' */hic_results/data/*/*merge*
grep 'rmdup' */hic_results/data/*/*merge*
grep 'cis_interaction' */hic_results/data/*/*merge*
grep 'cis_longRange' */hic_results/data/*/*merge*
grep 'trans_interaction' */hic_results/data/*/*merge*

grep 'Valid_interaction_pairs'$'\t' */hic_results/data/*/*.RSstat
grep 'Dangling_end_pairs' */hic_results/data/*/*.RSstat
grep 'Religation_pairs' */hic_results/data/*/*.RSstat
grep 'ycle' */hic_results/data/*/*.RSstat
