#!/bin/bash
# Submit process_sample.sh tasks first..last, keeping at most MAX of them running or queued.
#   submit_samples.sh 1 33
MAX=5
cd /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116
for i in $(seq $1 $2); do
    while [ $(qstat -u nergaon | grep -c process_sa || true) -ge $MAX ]; do sleep 60; done
    qsub -t $i process_sample.sh
done
