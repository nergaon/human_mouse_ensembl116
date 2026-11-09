#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -l hostname=sge1051
#$ -N submitter
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# Runs submit_samples.sh as a cluster job, so it does not depend on an interactive session.
#   qsub submitter_job.sh <first row> <last row>
cd /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116
./submit_samples.sh $1 $2
echo "submitted rows $1-$2"
