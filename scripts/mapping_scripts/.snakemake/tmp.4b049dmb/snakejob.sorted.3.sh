#!/bin/sh
# properties = {"type": "single", "rule": "sorted", "local": false, "input": ["/workdir/users/agk85/CDC2/mapping/B357-2/B357-2.bam"], "output": ["/workdir/users/agk85/CDC2/mapping/B357-2/B357-2.sorted.bam", "/workdir/users/agk85/CDC2/mapping/B357-2/done"], "wildcards": {"NAME": "B357-2"}, "params": {"n": "B357-2_sort", "NAME": "B357-2"}, "log": [], "threads": 1, "resources": {"mem_mb": 50}, "jobid": 3, "cluster": {}}
cd /local/workdir/users/agk85/CDC2/scripts/mapping && \
/usr/bin/python3.6 \
-m snakemake /workdir/users/agk85/CDC2/mapping/B357-2/done --snakefile /workdir/users/agk85/CDC2/scripts/mapping/snake_bwa_scaffold \
--force -j --keep-target-files --keep-remote \
--wait-for-files /local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb /workdir/users/agk85/CDC2/mapping/B357-2/B357-2.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules sorted --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/3.jobfinished" || (touch "/local/workdir/users/agk85/CDC2/scripts/mapping/.snakemake/tmp.4b049dmb/3.jobfailed"; exit 1)

