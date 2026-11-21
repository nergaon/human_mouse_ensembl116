#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -l hostname=sge1051
#$ -N rerun_lca6
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# Classification with the LeafCutter six-group test instead of the scipy ANOVA (user, 2026-11-21): LeafCutter for the
# 4 orthology variants x 5 cell types (run_lcanova6_116.sh array), then classification, downstream analyses, figures,
# tables and documents. Events, effects and per-pair LeafCutter runs are unchanged (rerun_liftover_6ds.sh).
#   qsub six/rerun_lcanova_6ds.sh
# Progress: ensembl116/logs/rerun_lcanova6.status; step logs rerunLC_<step>.log. Previous state: ensembl116/pre_lcanova_20261121/
set -u
BASE=/gpfs0/tals/projects/Analysis/human_mouse_exons
E=$BASE/ensembl116
S=$BASE/scripts/ensembl116
PY=/gpfs0/tals/projects/software/anaconda3/envs/spyder-env/bin/python
VPY=$S/.venv/bin/python
STATUS=$E/logs/rerun_lcanova6.status
: > $STATUS
step() { local name=$1; shift; echo "$(date '+%F %T') start $name" >> $STATUS
         "$@" > $E/logs/rerunLC_$name.log 2>&1; echo "$(date '+%F %T') exit $? $name" >> $STATUS; }
wait_jobs() { while qstat -u nergaon | awk 'NR>2{print $3}' | grep -q "^$1"; do sleep 60; done; }
missing() { local i=0; for v in id50 id70 id80 id50_rbh; do for c in CD4T CD8T NveB NK Mono; do i=$((i+1))
            f=$E/events/${v}_6ds_anova/leafcutter_anova6/$c.tsv
            [ -s $f ] && [ $(wc -l < $f) -gt 1000 ] || echo $i; done; done; }

cd $S/six
for v in id50 id70 id80 id50_rbh; do rm -rf $E/events/${v}_6ds_anova/leafcutter_anova6; done
# 1. LeafCutter six-group test (20 tasks); resubmit missing tasks once
qsub -N lca6 -t 1-20 run_lcanova6_116.sh >> $STATUS
sleep 30; wait_jobs lca6
for t in $(missing); do qsub -N lca6r -t $t run_lcanova6_116.sh >> $STATUS; done
sleep 30; wait_jobs lca6r
echo "LeafCutter six-group outputs missing: $(missing | wc -l)" >> $STATUS

# 2. classification and downstream analyses
for v in id50 id70 id80 id50_rbh; do step classify_$v $VPY classify_anova_116.py ${v}_6ds; done
step minor_usage $VPY minor_usage_116.py id50_6ds_anova
step single_cell $VPY single_cell_breakdown_116.py
step nagnag $VPY nagnag_test_116.py
step revision $VPY revision_analyses_116.py id50_6ds_anova
step conserved_options $VPY conserved_options_116.py
step facts $VPY facts_116.py

# 3. figures, tables, documents (the numbers in the text are updated afterwards)
step gene_panels bash make_gene_panels_6ds.sh
step figures_main $VPY figures_main_116.py
step figure_s3 $PY figure_classified_dpsi_116.py
step confounder_outputs $VPY confounder_outputs_116.py
step fig_junctions $VPY figure_status_by_junctions_116.py
step fig_minor $VPY figure_minor_usage_dpsi_116.py
step fig_cooc $VPY figure_status_cooccurrence_116.py
step fig_samplecorr $VPY figure_sample_correlation_116.py
step build_supp $PY build_supp_tables_116.py
step tables_b $VPY make_tables_b_116.py
step supp_docx $VPY make_supplementary_docx_116.py
step manuscript $VPY make_revised_docx_116.py $E/ref/Human_Mouse_Splicing_all_genes_HN44.docx $E/ref/hn44_unpacked $E/revision_6ds/manuscript
step letter $VPY response_letter_116.py
step ts2 $VPY ts2_responses_116.py
echo "$(date '+%F %T') done" >> $STATUS
