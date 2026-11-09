#!/usr/bin/env python3
"""
Species separation (sensitivity analysis, not in the manuscript): an event is species-separated in a cell type
when, for at least one of its junctions, the PSI of every human dataset (H1, H2, H3) differs from the PSI of every
mouse dataset (M1, M2) by at least a gap g (g = 0.1, 0.2): gap = min(H) - max(M) or min(M) - max(H).
Dataset PSI = mean of the per-sample PSI over samples with >= 20 cluster reads (as effect_size_analysis.py);
testable = PSI in all five datasets. Null: the same count for the other nine 3-vs-2 splits of the five datasets
(splits that mix the species); expected false = their mean, FDR = expected / observed.
Equivalent: every junction has PSI within a range of 0.1 across all five datasets.
Output: ensembl116/sensitivity/Table1_separation.docx / .tsv, species_separation_events.tsv
"""
import itertools
import os

import numpy as np
import pandas as pd

from table1_noncomplex_116 import to_docx, to_tsv

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50/'
OUT = B + 'sensitivity/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
DATASETS = ['H1', 'H2', 'H3', 'M1', 'M2']
HUMAN = ('H1', 'H2', 'H3')
GAPS = (0.1, 0.2)
MIN_READS = 20


def dataset_psi(counts, design, cell):
    cluster = counts.index.str.split(':').str[3]
    psi = {}
    for ds in DATASETS:
        s = design[(design['cell'] == cell) & (design['dataset'] == ds)].index.intersection(counts.columns)
        tot = counts[s].groupby(cluster).transform('sum')
        psi[ds] = (counts[s] / tot.where(tot >= MIN_READS)).mean(axis=1)
    psi = pd.DataFrame(psi)
    # testable: every junction of the cluster has a PSI in every dataset
    ok = psi.notna().all(axis=1).groupby(cluster).transform('all')
    psi = psi[ok.to_numpy()]
    return psi, psi.index.str.split(':').str[3]


def gap(psi, group):
    other = [d for d in DATASETS if d not in group]
    a, b = psi[list(group)], psi[other]
    return np.maximum(a.min(axis=1) - b.max(axis=1), b.min(axis=1) - a.max(axis=1))


def main():
    os.makedirs(OUT, exist_ok=True)
    counts = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(ED + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    status = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t').set_index(['cluster', 'cell'])
    ann = pd.read_csv(ED + 'classification/events_annotation.tsv', sep='\t').set_index('cluster')
    events, summary = [], {}
    for cell in CELLS:
        psi, clu = dataset_psi(counts, design, cell)
        splits = list(itertools.combinations(DATASETS, 3))
        g = pd.DataFrame({s: gap(psi, s).groupby(clu).max() for s in splits})
        rng = (psi.max(axis=1) - psi.min(axis=1)).groupby(clu).max()
        e = pd.DataFrame({'cell': cell, 'species_gap': g[HUMAN], 'max_range_all_datasets': rng})
        e['species_gap_rank'] = (g.ge(g[HUMAN], axis=0).sum(axis=1))   # 1 = no other split separates better
        events.append(e)
        s = {'testable': len(g), 'equivalent': int((rng <= 0.1).sum())}
        for t in GAPS:
            obs = int((g[HUMAN] >= t).sum())
            null = (g.drop(columns=[HUMAN]) >= t).sum()
            s[t] = dict(obs=obs, exp=null.mean(), max=int(null.max()))
            sep = g.index[g[HUMAN] >= t]
            st = status.reindex([(c, cell) for c in sep])
            s[t]['ds_hl'] = int((st['status_high_low'] == 'divergent').sum())
            s[t]['ds_high'] = int((st['status_high'] == 'divergent').sum())
        summary[cell] = s
    ev = pd.concat(events).rename_axis('cluster').reset_index()
    ev = ev.join(ann[['symbol', 'event_type']], on='cluster')
    ev = ev.join(status[['status_confidence']], on=['cluster', 'cell'])
    ev.sort_values(['cell', 'species_gap'], ascending=[True, False]).round(4) \
        .to_csv(OUT + 'species_separation_events.tsv', sep='\t', index=False)

    def row(label, f):
        return ('row', [label] + [f(summary[c]) for c in CELLS] + ['–'])

    def pct(n, s):
        return f'{n} ({100 * n / s["testable"]:.1f}%)'
    R = [('header', ['Comparison category', 'CD4T', 'CD8T', 'B', 'NK', 'Mono', 'Fibroblasts']),
         ('section', ['Human-Mouse alternative splicing events, species separation (PSI of each of the 3 human and '
                      '2 mouse datasets; fibroblasts: 1 dataset per species, not applicable)']),
         row('Testable events (PSI in all five datasets)', lambda s: str(s['testable']))]
    for t in GAPS:
        R += [row(f'Species-separated events, gap ≥ {t}', lambda s, t=t: pct(s[t]['obs'], s)),
              row(f'   expected by chance, mean [max] of the 9 other 3-vs-2 dataset splits',
                  lambda s, t=t: f"{s[t]['exp']:.0f} [{s[t]['max']}]"),
              row('   estimated false discovery rate, %',
                  lambda s, t=t: f"{100 * s[t]['exp'] / s[t]['obs']:.0f}" if s[t]['obs'] else '–'),
              row('   differentially spliced, high + low confidence (classification)',
                  lambda s, t=t: f"{s[t]['ds_hl']} ({100 * s[t]['ds_hl'] / max(s[t]['obs'], 1):.0f}%)"),
              row('   differentially spliced, high confidence (classification)',
                  lambda s, t=t: f"{s[t]['ds_high']} ({100 * s[t]['ds_high'] / max(s[t]['obs'], 1):.0f}%)")]
    R.append(row('Events with PSI range ≤ 0.1 across all five datasets', lambda s: pct(s['equivalent'], s)))
    to_tsv(R, OUT + 'Table1_separation.tsv')
    to_docx(R, 'Table 1, species separation. An event is species-separated in a cell type when, for at least one '
               'junction, the PSI of every human dataset differs from the PSI of every mouse dataset by at least the '
               'gap. The expected number by chance is the same count for the nine other ways of splitting the five '
               'datasets into groups of three and two, which mix the species. Percentages in the classification rows '
               'are of the species-separated events.', OUT + 'Table1_separation.docx')
    print(pd.read_csv(OUT + 'Table1_separation.tsv', sep='\t', header=None).to_string(index=False, header=False))


if __name__ == '__main__':
    main()
