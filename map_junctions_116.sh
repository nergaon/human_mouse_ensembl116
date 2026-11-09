#!/bin/bash
# Orthologous junctions of every dataset, for each Ensembl 116 orthology variant.
# Human tables are filtered to junctions whose two ends are orthologous boundaries (position_h);
# mouse tables are remapped from position_m to the orthologous human coordinates.
#   map_junctions_116.sh            (all variants)
#   map_junctions_116.sh id50       (one variant)
# Output: ensembl116/junctions/<variant>/<study>_orthologs_junctions.tsv
set -euo pipefail
BASE=/gpfs0/tals/projects/Analysis/human_mouse_exons
SCRIPT=$BASE/scripts/ensembl116/filter_ortholog_junctions.py
VARIANTS=${1:-"id50 id70 id80 id50_rbh"}

declare -A HUMAN=(
  [GSE60424]=$BASE/GSE60424/leafcutter_0.2.9/GSE60424_JSR_junction_counts_no_strand.tsv
  [GSE115736]=$BASE/GSE115736/leafcutter_0.2.9/GSE115736_JSR_junction_counts_no_strand.tsv
  [GSE118165]=$BASE/GSE118165/leafcutter_0.2.9/GSE118165_JSR_junction_counts_no_strand.tsv
  [GSE107011]=$BASE/GSE107011/leafcutter_0.2.9/GSE107011_JSR_junction_counts_no_strand.tsv
  [GSE133383]=$BASE/GSE133383/leafcutter_0.2.9/GSE133383_JSR_junction_counts_no_strand.tsv  # tissue control
  [GSE143213]=$BASE/GSE143213/leafcutter_0.2.9/GSE143213_JSR_junction_counts_no_strand.tsv  # age control
  [EMTAB5919H]=$BASE/fibroblast/E-MTAB-5919/leafcutter_0.2.9/E-MTAB-5919_H_JSR_junction_counts_no_strand.tsv
)
declare -A MOUSE=(
  [GSE116177]=$BASE/GSE116177/leafcutter_0.2.9/GSE116177_JSR_junction_counts_no_strand.tsv
  [GSE287527]=$BASE/GSE287527/leafcutter_0.2.9/GSE287527_JSR_junction_counts_no_strand.tsv
  [GSE261448]=$BASE/GSE261448/leafcutter_0.2.9/GSE261448_JSR_junction_counts_no_strand.tsv
  [ENCODE4]=$BASE/ENCODE4/leafcutter_0.2.9/ENCODE4_JSR_junction_counts_no_strand.tsv
  [GSE313371]=$BASE/GSE313371/leafcutter_0.2.9/GSE313371_JSR_junction_counts_no_strand.tsv
  [GSE182858]=$BASE/GSE182858/leafcutter_0.2.9/GSE182858_JSR_junction_counts_no_strand.tsv
  [GSE124829]=$BASE/GSE124829/leafcutter_0.2.9/GSE124829_JSR_junction_counts_no_strand.tsv  # sex / age controls
  [ENCODEB]=$BASE/ENCODEB/leafcutter_0.2.9/ENCODEB_JSR_junction_counts_no_strand.tsv  # sex control
  [GSE233347]=$BASE/GSE233347/leafcutter_0.2.9/GSE233347_JSR_junction_counts_no_strand.tsv  # tissue control
  [GSE184496]=$BASE/GSE184496/leafcutter_0.2.9/GSE184496_JSR_junction_counts_no_strand.tsv  # age control
  [EMTAB5919M]=$BASE/fibroblast/E-MTAB-5919/leafcutter_0.2.9/E-MTAB-5919_M_JSR_junction_counts_no_strand.tsv
)

for v in $VARIANTS; do
    out=$BASE/ensembl116/junctions/$v
    points=$BASE/ensembl116/orthology/unique_points_116_$v.txt
    mkdir -p $out
    for study in "${!HUMAN[@]}"; do
        [[ -s ${HUMAN[$study]} ]] || { echo "missing ${HUMAN[$study]}"; continue; }
        python3 $SCRIPT --points $points --mode filter-only --source-column position_h \
            --inputs ${HUMAN[$study]} --outputs $out/${study}_orthologs_junctions.tsv
    done
    for study in "${!MOUSE[@]}"; do
        [[ -s ${MOUSE[$study]} ]] || { echo "missing ${MOUSE[$study]}"; continue; }
        python3 $SCRIPT --points $points --mode remap --source-column position_m --target-column position_h \
            --inputs ${MOUSE[$study]} --outputs $out/${study}_orthologs_junctions.tsv
    done
done
