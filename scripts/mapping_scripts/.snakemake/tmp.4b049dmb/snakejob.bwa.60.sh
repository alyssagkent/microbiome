#!/bin/sh
# properties = {"type": "single", "rule": "bwa", "local": false, "input": ["/workdir/data/CDC/metagenomes/merged/unzip/B370-4.1.fastq", "/workdir/data/CDC/metagenomes/merged/unzip/B370-4.2.fastq"], "output": ["/workdir/users/agk85/CDC2/mapping/B370-4/B370-4.sam"], "wildcards": {"NAME": "B370-4"}, "params": {"REF": "/workdir/users/agk85/CDC2/hic/mapping/references/B370-4_scaffold", "MAPFOLDER": "/workdir/users/agk85/CDC2/mapping/B370-4", "n": "B370-4_bwa"}, "log": [], "threads": 1, "resources": {"mem_mb": 40}, "jobid": 60, "cluster": {}}
cd /local/workdir/users/agk85/CDC2/scripts/mapping && \
/usr/bin/python3.6 \
-m snakemake /workdir/users/agk85/CDC2/mapping/B370-4/B370-4.sam --snakefile /workdir/users/agk85/CDC2/scripts/mapping/snake_bwa_scaffold \
--force -j --keep-target-files --keep-remote \
--wait-for-files /local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb /workdir/data/CDC/metagenomes/merged/unzip/B370-4.1.fastq /workdir/data/CDC/metagenomes/merged/unzip/B370-4.2.fastq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules bwa --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/60.jobfinished" || (touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/60.jobfailed"; exit 1)

