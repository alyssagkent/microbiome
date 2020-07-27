#!/bin/sh
# properties = {"type": "single", "rule": "bwa", "local": false, "input": ["/workdir/data/CDC/metagenomes/merged/unzip/US3-16.1.fastq", "/workdir/data/CDC/metagenomes/merged/unzip/US3-16.2.fastq"], "output": ["/workdir/users/agk85/CDC2/mapping/US3-16/US3-16.sam"], "wildcards": {"NAME": "US3-16"}, "params": {"REF": "/workdir/users/agk85/CDC2/hic/mapping/references/US3-16_scaffold", "MAPFOLDER": "/workdir/users/agk85/CDC2/mapping/US3-16", "n": "US3-16_bwa"}, "log": [], "threads": 1, "resources": {"mem_mb": 40}, "jobid": 62, "cluster": {}}
cd /local/workdir/users/agk85/CDC2/scripts/mapping && \
/usr/bin/python3.6 \
-m snakemake /workdir/users/agk85/CDC2/mapping/US3-16/US3-16.sam --snakefile /workdir/users/agk85/CDC2/scripts/mapping/snake_bwa_scaffold \
--force -j --keep-target-files --keep-remote \
--wait-for-files /local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb /workdir/data/CDC/metagenomes/merged/unzip/US3-16.1.fastq /workdir/data/CDC/metagenomes/merged/unzip/US3-16.2.fastq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules bwa --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/62.jobfinished" || (touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/62.jobfailed"; exit 1)

