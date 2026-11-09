#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -pe shared 4
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# LeafCutter differential splicing, one task per line of jobs.txt (from leafcutter_design_116.py).
# Same parameters as leafcutter_0.2.9_HN5.sh: -i = -g = max(2, floor(0.6 * smaller group)),
# min_coverage default (20). Cluster gene labels from the Ensembl 116 exon file.
#   qsub -t 1-<n jobs> run_leafcutter_116.sh <leafcutter dir containing jobs.txt> <AS_clusters_value.txt>
# Use the env without 'conda activate': its r-base activation hook runs 'R CMD javareconf', which
# rewrites lib/R/etc/ldpaths and Makeconf; parallel tasks racing on that left both files empty.
LEAFCUTTER_ENV=/gpfs0/tals/projects/software/anaconda3/envs/leafcutter
export PATH=$LEAFCUTTER_ENV/bin:$PATH
LEAFCUTTER=/gpfs0/tals/projects/software/leafcutter_0.2.9/leafcutter/scripts/leafcutter_ds.R
EXON_FILE=/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/ref/ensembl116_human_exons_leafcutter.txt.gz

folder=$(sed -n "${SGE_TASK_ID}p" $1/jobs.txt)
counts=$2
cd $folder
min_group_size=$(cut -f2 groups_file.txt | sort | uniq -c | awk '{print $1}' | sort -n | head -n1)
threshold=$(( min_group_size * 60 / 100 ))
[ "$threshold" -lt 2 ] && threshold=2
echo "$folder min group $min_group_size threshold $threshold"
Rscript $LEAFCUTTER --num_threads 4 -i $threshold -g $threshold --exon_file=$EXON_FILE $counts groups_file.txt
