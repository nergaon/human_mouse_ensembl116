#!/usr/bin/env python3
"""
Calibration of a pooled LeafCutter test (not part of the manuscript).

For each immune cell type and each of the ten 3-vs-2 splits of the five datasets (H1, H2, H3, M1, M2), all samples
of the three datasets of one group are compared with all samples of the two datasets of the other group in one
LeafCutter run. The human-mouse split (H1 H2 H3 vs M1 M2) is the pooled species test; the nine mixed splits show how
many events such a pooled test calls significant when the groups mix the species (dataset differences only).
Samples: the 'original' sample sets of the main analysis (leafcutter/sample_design.tsv of id50).

Output: ensembl116/events/id50_pooled/leafcutter/<split>/<cell>/groups_file.txt and leafcutter/jobs.txt
Then:   qsub -N lc_pooled -t 1-50 ../run_leafcutter_116.sh <out>/leafcutter <id50>/AS_clusters_value.txt
"""
import itertools
import os

import pandas as pd

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/'
OUT = BASE + 'id50_pooled/leafcutter/'
DATASETS = ['H1', 'H2', 'H3', 'M1', 'M2']
CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']


def main():
    design = pd.read_csv(BASE + 'id50/leafcutter/sample_design.tsv', sep='\t')
    design = design[design['dataset'].isin(DATASETS)]
    jobs = []
    for a in itertools.combinations(DATASETS, 3):
        b = tuple(d for d in DATASETS if d not in a)
        name = ''.join(a) + '_vs_' + ''.join(b)
        for cell in CELLS:
            d = design[design['cell'] == cell]
            ga, gb = d[d['dataset'].isin(a)]['sample'], d[d['dataset'].isin(b)]['sample']
            path = OUT + f'{name}/{cell}'
            os.makedirs(path, exist_ok=True)
            pd.concat([pd.DataFrame({'s': ga, 'g': 'A'}), pd.DataFrame({'s': gb, 'g': 'B'})]) \
                .to_csv(path + '/groups_file.txt', sep='\t', header=False, index=False)
            jobs.append(path)
            print(name, cell, len(ga), len(gb))
    with open(OUT + 'jobs.txt', 'w') as f:
        f.write('\n'.join(jobs) + '\n')
    print(len(jobs), 'runs')


if __name__ == '__main__':
    main()
