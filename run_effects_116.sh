#!/bin/bash
# Effect sizes / equivalence classes and their summary for every LeafCutter comparison of one variant.
#   run_effects_116.sh id50
# Output: ensembl116/events/<variant>/effects/
set -euo pipefail
E=/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/$1
cd /gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116
args=()
for d in $E/leafcutter/original/*/; do args+=(--comparison "$(basename $d)=$d"); done
for d in $E/leafcutter/matched/*/; do args+=(--comparison "matched_$(basename $d)=$d"); done
python3 effect_size_analysis.py --counts $E/AS_clusters_value.txt "${args[@]}" --out-dir $E/effects > $E/effects.log 2>&1
python3 summarize_comparisons_116.py $E/effects
