cd /workdir/users/agk85/CDC2/das

cat B314*/*_bintable.txt | cut -f5 | sort | uniq > B314_taxa.txt
cat B316*/*_bintable.txt | cut -f5 | sort | uniq > B316_taxa.txt
cat B320*/*_bintable.txt | cut -f5 | sort | uniq > B320_taxa.txt
cat B331*/*_bintable.txt | cut -f5 | sort | uniq > B331_taxa.txt
cat B335*/*_bintable.txt | cut -f5 | sort | uniq > B335_taxa.txt
cat B357*/*_bintable.txt | cut -f5 | sort | uniq > B357_taxa.txt
cat B370*/*_bintable.txt | cut -f5 | sort | uniq > B370_taxa.txt
cat US3*/*_bintable.txt | cut -f5 | sort | uniq > US3_taxa.txt
cat US8*/*_bintable.txt | cut -f5 | sort | uniq > US8_taxa.txt
cat */*_bintable.txt | cut -f5 | sort | uniq > all_taxa.txt
