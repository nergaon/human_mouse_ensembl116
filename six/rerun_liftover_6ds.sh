#!/bin/bash
#$ -cwd
#$ -V
#$ -q tals.q
#$ -l hostname=sge1051
#$ -N rerun_lo6
#$ -o /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
#$ -e /gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/logs/
# Six-dataset rerun after the liftOver correction of the orthologous boundary pairs (liftover/liftover_boundaries.py,
# corrected pairs installed as orthology/unique_points_116_<variant>.txt) and with the splice-site QC
# (qc_splice_sites_116.py). One cluster job that submits and waits for the LeafCutter arrays (sge1051 is a submit host).
# Previous outputs: ensembl116/pre_liftover_20261118/.
#   qsub six/rerun_liftover_6ds.sh
# Progress: ensembl116/logs/rerun_liftover6.status (one line per step with its exit code); step logs rerun6_<step>.log
set -u
BASE=/gpfs0/tals/projects/Analysis/human_mouse_exons
E=$BASE/ensembl116
S=$BASE/scripts/ensembl116
PY=/gpfs0/tals/projects/software/anaconda3/envs/spyder-env/bin/python   # networkx
VPY=$S/.venv/bin/python
STATUS=$E/logs/rerun_liftover6.status
VARIANTS="id50 id70 id80 id50_rbh"
E6=$E/events/id50_6ds
A6=$E/events/id50_6ds_anova
: > $STATUS

step() {  # step <name> <command...>
    local name=$1; shift
    echo "$(date '+%F %T') start $name" >> $STATUS
    "$@" > $E/logs/rerun6_$name.log 2>&1
    echo "$(date '+%F %T') exit $? $name" >> $STATUS
}
wait_jobs() { while qstat -u nergaon | awk 'NR>2{print $3}' | grep -q "^$1"; do sleep 60; done; }
lc_failed() { local i=0; while read d; do i=$((i+1)); [ -s $d/leafcutter_ds_cluster_significance.txt ] || echo $i; done < $1; }

# 1. orthologous junction tables of all studies (the m3 studies are remapped by build_events_6ds.py)
cd $S
rm -f $E/m3/junctions/*_orthologs_junctions.tsv
rm -rf $E/m3/junctions/id70 $E/m3/junctions/id80 $E/m3/junctions/id50_rbh
step map_junctions bash $S/map_junctions_116.sh

# 2. events, splice-site QC
cd $S/six
for v in $VARIANTS; do
    rm -rf $E/events/${v}_6ds
    step build_$v $VPY build_events_6ds.py $v
    step qc_$v $VPY qc_splice_sites_116.py $v
done

# 3. LeafCutter: 15 dataset pairs + fibroblasts; sex / age / tissue controls
L=$E6/leafcutter
qsub -N lc6_main -t 1-$(wc -l < $L/jobs.txt) $S/run_leafcutter_116.sh $L $E6/AS_clusters_value.txt >> $STATUS
J=$E/junctions/id50
step conf_counts $PY confounders_116.py counts GSE107011:$J/GSE107011_orthologs_junctions.tsv \
    GSE133383:$J/GSE133383_orthologs_junctions.tsv GSE143213:$J/GSE143213_orthologs_junctions.tsv \
    ENCODEB:$J/ENCODEB_orthologs_junctions.tsv GSE233347:$J/GSE233347_orthologs_junctions.tsv \
    GSE184496:$J/GSE184496_orthologs_junctions.tsv GSE116177:$J/GSE116177_orthologs_junctions.tsv \
    PRJNA797236:$J/PRJNA797236_orthologs_junctions.tsv
step conf_design $PY confounders_116.py design
C=$A6/leafcutter/confounders
qsub -N lc6_conf -t 1-$(wc -l < $C/jobs.txt) $S/run_leafcutter_116.sh $C/ $A6/AS_clusters_value_confounders.txt >> $STATUS
wait_jobs lc6_
for t in $(lc_failed $L/jobs.txt); do qsub -N lc6r -t $t $S/run_leafcutter_116.sh $L $E6/AS_clusters_value.txt >> $STATUS; done
for t in $(lc_failed $C/jobs.txt); do qsub -N lc6r -t $t $S/run_leafcutter_116.sh $C/ $A6/AS_clusters_value_confounders.txt >> $STATUS; done
wait_jobs lc6r
echo "LeafCutter missing outputs: main $(lc_failed $L/jobs.txt | wc -l), controls $(lc_failed $C/jobs.txt | wc -l)" >> $STATUS

# 4. effect sizes, per-comparison status, classification
cd $S
args=(); for d in $L/original/*/; do args+=(--comparison "$(basename $d)=$d"); done
step effects bash -c "python3 effect_size_analysis.py --counts $E6/AS_clusters_value.txt ${args[*]} --out-dir $E6/effects && python3 summarize_comparisons_116.py $E6/effects"
cargs=(); for d in $C/*/; do cargs+=(--comparison "$(basename $d)=$d"); done
step effects_conf python3 effect_size_analysis.py --counts $A6/AS_clusters_value_confounders.txt "${cargs[@]}" --out-dir $E6/effects_confounders
step classify_events $PY classify_events_116.py id50_6ds
cd $S/six
step single_status $VPY single_status_116.py
for v in $VARIANTS; do step classify_$v $VPY classify_anova_116.py ${v}_6ds; done

# 5. downstream analyses (figures, tables and text are rebuilt interactively afterwards)
step minor_usage $VPY minor_usage_116.py id50_6ds_anova
step single_cell $VPY single_cell_breakdown_116.py
step nagnag $VPY nagnag_test_116.py
step conf_summary $VPY confounders_116.py summary
step revision $VPY revision_analyses_116.py id50_6ds_anova
step facts $VPY facts_116.py
echo "$(date '+%F %T') done" >> $STATUS
