#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -pe shared 8
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# One task per (orthology variant, cell type): LeafCutter six-group test (leafcutter_anova6_116.R).
#   qsub -t 1-20 run_lcanova6_116.sh
# The leafcutter env is used without 'conda activate' (its activation hook races on lib/R/etc in parallel tasks).
export PATH=/gpfs0/tals/projects/software/anaconda3/envs/leafcutter/bin:$PATH
VARIANTS=(id50 id70 id80 id50_rbh)
CELLS=(CD4T CD8T NveB NK Mono)
i=$((SGE_TASK_ID - 1))
v=${VARIANTS[$((i / 5))]}
c=${CELLS[$((i % 5))]}
cd /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116/six
Rscript leafcutter_anova6_116.R $v $c 8
