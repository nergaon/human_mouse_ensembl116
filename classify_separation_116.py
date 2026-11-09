#!/usr/bin/env python3
"""
Human-mouse event classification by dataset separation (alternative to classify_events_116.py; used by the
revision in ensembl116/revision_separation/). Per cell type, from the dataset PSI (mean of the per-sample PSI over
samples with >= 20 cluster reads) of the three human (H1-H3) and two mouse (M1, M2) datasets:
  testable   every junction of the event has a PSI in all five datasets
  gap        max over junctions of the distance between the human and the mouse PSI ranges
             (min(H) - max(M) or min(M) - max(H); > 0 only if every human dataset is on one side of every mouse one)
  range      max over junctions of the PSI range over all five datasets (max - min)
  divergent (DS)   high: gap >= 0.2   high + low: gap >= 0.1
  conserved (equivalent)   high: range <= 0.05   high + low: range <= 0.1
  inconclusive     all other testable events
Null for the divergent class: the same gap for the nine other 3-vs-2 splits of the five datasets (columns
gap_rank, null_splits_gap_*). The per-comparison columns (cross_*, within_*, status_<pair>, DS_within_species) are kept
from classify_events_116.py (LeafCutter + effect sizes) for the analyses that use them.

Usage: classify_separation_116.py <variant>      (e.g. id50)
Output: ensembl116/events/<variant>_sep/  (classification/, revision/, and links to the data of <variant>)
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
DATASETS = ['H1', 'H2', 'H3', 'M1', 'M2']
HUMAN = ('H1', 'H2', 'H3')
GAP_HIGH, GAP_LOW = 0.2, 0.1
RANGE_HIGH, RANGE_LOW = 0.05, 0.1
MIN_READS = 20
LINKS = ['AS_clusters_value.txt', 'AS_clusters_value_confounders.txt', 'junction_presence.tsv', 'effects',
         'effects_confounders', 'leafcutter']
LEVELS = ('high', 'high_low')


def dataset_psi(counts, design, cell):
    cluster = counts.index.str.split(':').str[3]
    psi = {}
    for ds in DATASETS:
        s = design[(design['cell'] == cell) & (design['dataset'] == ds)].index.intersection(counts.columns)
        tot = counts[s].groupby(cluster).transform('sum')
        psi[ds] = (counts[s] / tot.where(tot >= MIN_READS)).mean(axis=1)
    psi = pd.DataFrame(psi)
    ok = psi.notna().all(axis=1).groupby(cluster).transform('all')
    psi = psi[ok.to_numpy()]
    return psi, psi.index.str.split(':').str[3]


def gap(psi, group):
    other = [d for d in DATASETS if d not in group]
    a, b = psi[list(group)], psi[other]
    return np.maximum(a.min(axis=1) - b.max(axis=1), b.min(axis=1) - a.max(axis=1))


def separation(variant):
    ed = BASE + f'events/{variant}/'
    counts = pd.read_csv(ed + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(ed + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    out = []
    for cell in CELLS:
        psi, clu = dataset_psi(counts, design, cell)
        g = pd.DataFrame({s: gap(psi, s).groupby(clu).max() for s in itertools.combinations(DATASETS, 3)})
        null = g.drop(columns=[HUMAN])
        out.append(pd.DataFrame({
            'cell': cell, 'species_gap': g[HUMAN],
            'psi_range': (psi.max(axis=1) - psi.min(axis=1)).groupby(clu).max(),
            'gap_rank': g.ge(g[HUMAN], axis=0).sum(axis=1),
            'null_splits_gap_low': null.ge(GAP_LOW).sum(axis=1),       # of the 9 other splits
            'null_splits_gap_high': null.ge(GAP_HIGH).sum(axis=1)}).rename_axis('cluster').reset_index())
    return pd.concat(out, ignore_index=True)


def main():
    variant = sys.argv[1]
    new = BASE + f'events/{variant}_sep/'
    os.makedirs(new + 'classification', exist_ok=True)
    os.makedirs(new + 'revision', exist_ok=True)
    for name in LINKS:
        if os.path.exists(BASE + f'events/{variant}/{name}') and not os.path.lexists(new + name):
            os.symlink(BASE + f'events/{variant}/{name}', new + name)
    old = BASE + f'events/{variant}/classification/'
    ann = pd.read_csv(old + 'events_annotation.tsv', sep='\t')
    ann.to_csv(new + 'classification/events_annotation.tsv', sep='\t', index=False)
    status = pd.read_csv(old + 'event_status_by_cell.tsv', sep='\t')
    sep = separation(variant)
    status = status.drop(columns=['status_confidence', 'status_high', 'status_high_low']) \
        .merge(sep, on=['cell', 'cluster'], how='left')
    g, r, t = status['species_gap'], status['psi_range'], status['species_gap'].notna()
    status['status_confidence'] = np.select(
        [~t, g >= GAP_HIGH, g >= GAP_LOW, r <= RANGE_HIGH, r <= RANGE_LOW],
        ['not testable', 'divergent high', 'divergent low', 'conserved high', 'conserved low'], 'inconclusive')
    status['status_high'] = np.select([~t, g >= GAP_HIGH, r <= RANGE_HIGH],
                                      ['not testable', 'divergent', 'conserved'], 'inconclusive')
    status['status_high_low'] = np.select([~t, g >= GAP_LOW, r <= RANGE_LOW],
                                          ['not testable', 'divergent', 'conserved'], 'inconclusive')
    status.to_csv(new + 'classification/event_status_by_cell.tsv', sep='\t', index=False)
    summary = pd.concat([status.groupby('cell')[f'status_{n}'].value_counts().unstack(fill_value=0)
                         .add_prefix(f'{n}_') for n in LEVELS], axis=1)
    for n in ('divergent', 'conserved'):
        summary[f'high_low_{n}_DS_within_species'] = (status['DS_within_species'] &
                                                     (status['status_high_low'] == n)).groupby(status['cell']).sum()
    # expected false divergent events: mean count over the nine splits that mix the species
    for lab, col in (('high', 'null_splits_gap_high'), ('high_low', 'null_splits_gap_low')):
        summary[f'{lab}_divergent_expected'] = (status[col] / 9).groupby(status['cell']).sum().round(1)
    summary.to_csv(new + 'classification/summary_by_cell.tsv', sep='\t')
    print(summary.to_string())


if __name__ == '__main__':
    main()
