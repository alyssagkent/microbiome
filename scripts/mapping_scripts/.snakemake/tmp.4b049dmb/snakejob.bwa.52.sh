#!/bin/sh
# properties = {"type": "single", "rule": "bwa", "local": false, "input": ["/workdir/data/CDC/metagenomes/merged/unzip/US8-1.1.fastq", "/workdir/data/CDC/metagenomes/merged/unzip/US8-1.2.fastq"], "output": ["/workdir/users/agk85/CDC2/mapping/US8-1/US8-1.sam"], "wildcards": {"NAME": "US8-1"}, "params": {"REF": "/workdir/users/agk85/CDC2/hic/mapping/references/US8-1_scaffold", "MAPFOLDER": "/workdir/users/agk85/CDC2/mapping/US8-1", "n": "US8-1_bwa"}, "log": [], "threads": 1, "resources": {"mem_mb": 40}, "jobid": 52, "cluster": {}}
cd /local/workdir/users/agk85/CDC2/scripts/mapping && \
/usr/bin/python3.6 \
-m snakemake /workdir/users/agk85/CDC2/mapping/US8-1/US8-1.sam --snakefile /workdir/users/agk85/CDC2/scripts/mapping/snake_bwa_scaffold \
--force -j --keep-target-files --keep-remote \
--wait-for-files /local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb /workdir/data/CDC/metagenomes/merged/unzip/US8-1.1.fastq /workdir/data/CDC/metagenomes/merged/unzip/US8-1.2.fastq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules bwa --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/52.jobfinished" || (touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/52.jobfailed"; exit 1)

