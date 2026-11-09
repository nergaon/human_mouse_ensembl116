#!/bin/bash
# Downstream of map_junctions_116.sh for one orthology variant:
#   shared orthologous AS events -> LeafCutter design -> LeafCutter array job
#   run_downstream_116.sh id50
# Output: ensembl116/events/<variant>/ (AS_clusters_value.txt, leafcutter/<set>/<A>_<B>/<cell>/)
set -euo pipefail
v=$1
BASE=/gpfs0/tals/projects/Analysis/human_mouse_exons
J=$BASE/ensembl116/junctions/$v
OUT=$BASE/ensembl116/events/$v
PY=/gpfs0/tals/projects/software/anaconda3/envs/spyder-env/bin/python  # has networkx
cd $BASE/scripts/ensembl116
mkdir -p $OUT/leafcutter

# event universe from the immune datasets; fibroblast counts added for the same junctions
# (PBMC, GSE107011 vs GSE182858, dropped: mouse PBMC alignment 17-51%, mostly rRNA)
$PY shared_junctions_116.py \
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
    --out-dir $OUT 2>&1 | grep -v -E 'Processing chromosome|initial components'

python3 leafcutter_design_116.py --counts $OUT/AS_clusters_value.txt --out-dir $OUT/leafcutter
n=$(wc -l < $OUT/leafcutter/jobs.txt)
qsub -N lc_$v -t 1-$n run_leafcutter_116.sh $OUT/leafcutter $OUT/AS_clusters_value.txt
