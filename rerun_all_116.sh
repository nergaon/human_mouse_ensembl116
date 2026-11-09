#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -l hostname=sge1051
#$ -N rerun_all116
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# Full rerun downstream of orthology/unique_points_116_<variant>.txt (written by orthology_points_116.py),
# as one cluster job that submits and waits for the LeafCutter arrays (runs on sge1051, a submit host).
# Added 2026-11-02 after the boundary-pairing fix (mutual best) in orthology_points_116.py.
#   qsub rerun_all_116.sh
# Progress: ensembl116/logs/rerun_all116.status (one line per step with its exit code).
set -u
BASE=/gpfs0/tals/projects/Analysis/human_mouse_exons
E116=$BASE/ensembl116
S=$BASE/scripts/ensembl116
PY=/gpfs0/tals/projects/software/anaconda3/envs/spyder-env/bin/python   # networkx (shared junctions)
VPY=$S/.venv/bin/python                                                 # gseapy, docx, openpyxl
STATUS=$E116/logs/rerun_all116.status
VARIANTS="id50 id70 id80 id50_rbh"
cd $S
: > $STATUS

step() {  # step <name> <command...>: run, log exit code, keep going
    local name=$1; shift
    echo "$(date '+%F %T') start $name" >> $STATUS
    "$@" > $E116/logs/rerun_$name.log 2>&1
    echo "$(date '+%F %T') exit $? $name" >> $STATUS
}

wait_jobs() {  # wait until no job named $1* is queued or running
    while qstat -u nergaon | awk 'NR>2{print $3}' | grep -q "^$1"; do sleep 60; done
}

lc_failed() {  # tasks of a jobs.txt without LeafCutter output
    local n=0 i=0
    while read d; do i=$((i+1)); [ -s $d/leafcutter_ds_cluster_significance.txt ] || { n=$((n+1)); echo $i; }; done < $1
}

# 1. orthologous junction tables
step map_junctions bash map_junctions_116.sh

# 2. events and LeafCutter design per variant; LeafCutter arrays
for v in $VARIANTS; do
    OUT=$E116/events/$v; J=$E116/junctions/$v
    rm -rf $OUT; mkdir -p $OUT/leafcutter $OUT/revision
    step shared_$v $PY shared_junctions_116.py \
        --table H1:human:$J/GSE60424_orthologs_junctions.tsv \
        --table H2:human:$J/GSE115736_orthologs_junctions.tsv \
        --table H3:human:$J/GSE118165_orthologs_junctions.tsv \
        --table M1:mouse:$J/GSE116177_orthologs_junctions.tsv \
        --table M2:mouse:$J/GSE287527_orthologs_junctions.tsv \
        --table M2:mouse:$J/GSE261448_orthologs_junctions.tsv \
        --table M2:mouse:$J/ENCODE4_orthologs_junctions.tsv \
        --table M2:mouse:$J/GSE313371_orthologs_junctions.tsv \
        --extra FH:human:$J/EMTAB5919H_orthologs_junctions.tsv \
        --extra FM:mouse:$J/EMTAB5919M_orthologs_junctions.tsv \
        --out-dir $OUT
    step design_$v $PY leafcutter_design_116.py --counts $OUT/AS_clusters_value.txt --out-dir $OUT/leafcutter
    n=$(wc -l < $OUT/leafcutter/jobs.txt)
    qsub -N lc_$v -t 1-$n run_leafcutter_116.sh $OUT/leafcutter $OUT/AS_clusters_value.txt >> $STATUS
done

# 3. sex / age / tissue controls (id50 events)
J=$E116/junctions/id50
step conf_counts $PY confounders_116.py counts GSE107011:$J/GSE107011_orthologs_junctions.tsv \
    GSE133383:$J/GSE133383_orthologs_junctions.tsv GSE143213:$J/GSE143213_orthologs_junctions.tsv \
    ENCODEB:$J/ENCODEB_orthologs_junctions.tsv GSE233347:$J/GSE233347_orthologs_junctions.tsv \
    GSE184496:$J/GSE184496_orthologs_junctions.tsv
step conf_design $PY confounders_116.py design
C=$E116/events/id50/leafcutter/confounders
qsub -N lc_conf -t 1-$(wc -l < $C/jobs.txt) run_leafcutter_116.sh $C/ $E116/events/id50/AS_clusters_value_confounders.txt >> $STATUS

# 4. wait; resubmit failed LeafCutter tasks once
wait_jobs lc_
for v in $VARIANTS conf; do
    if [ $v = conf ]; then jobs=$C/jobs.txt; counts=$E116/events/id50/AS_clusters_value_confounders.txt; dir=$C/
    else jobs=$E116/events/$v/leafcutter/jobs.txt; counts=$E116/events/$v/AS_clusters_value.txt; dir=$E116/events/$v/leafcutter; fi
    for t in $(lc_failed $jobs); do qsub -N lc2_$v -t $t run_leafcutter_116.sh $dir $counts >> $STATUS; done
done
wait_jobs lc2_
for v in $VARIANTS; do echo "leafcutter $v failed tasks: $(lc_failed $E116/events/$v/leafcutter/jobs.txt | tr '\n' ' ')" >> $STATUS; done
echo "leafcutter conf failed tasks: $(lc_failed $C/jobs.txt | tr '\n' ' ')" >> $STATUS

# 5. effect sizes, classification, rank consistency per variant
for v in $VARIANTS; do
    step effects_$v bash run_effects_116.sh $v
    step classify_$v $PY classify_events_116.py $v
    step rank_$v $VPY rank_consistency_116.py $v
done
ED=$E116/events/id50
A=(); for d in $C/*/; do A+=(--comparison "$(basename $d)=$d"); done
step conf_effects python3 effect_size_analysis.py --counts $ED/AS_clusters_value_confounders.txt "${A[@]}" \
    --out-dir $ED/effects_confounders

# 6. id50 analyses (revision_analyses ends with the enrichment, which needs internet; rerun it interactively if
#    it fails - the other revision tables are written before it)
step minor_usage $VPY minor_usage_116.py id50
step single_cell $VPY single_cell_breakdown_116.py
step nagnag $VPY nagnag_test_116.py
step revision_analyses $VPY revision_analyses_116.py id50
step species_model $PY species_model_116.py id50
step conf_summary $VPY confounders_116.py summary
step facts $VPY facts_116.py

# 7. figures and supplementary tables (gene figures and the manuscript text are updated separately)
step figure1 $VPY figure1_116.py /gpfs0/tals/users/nergaon/home/.claude/uploads/local_8a84ed4e-ab88-45b4-b7d9-0727a9d85737/figure1_A4_v6.pdf
step figures_main $PY figures_main_116.py
step figure_dpsi_ci $VPY figure_dpsi_ci_116.py
step figure_classified $VPY figure_classified_dpsi_116.py
step confounder_outputs $VPY confounder_outputs_116.py
step supp_tables $VPY build_supp_tables_116.py
step supp_docx $VPY make_supplementary_docx_116.py
echo "$(date '+%F %T') ALL DONE" >> $STATUS
