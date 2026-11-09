#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -pe shared 32
#$ -N exon_orth116
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# qsub run_exon_orthology_116.sh
set -e
cd /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116
python3 exon_orthology_116.py --threads 32
python3 orthology_points_116.py
