#!/usr/bin/env python3
"""
Hierarchical first step (as in Ner-Gaon et al. 2023, Front Immunol 14:1116392: one test across all groups first,
pairwise tests only for the events significant in it): one-way ANOVA of junction PSI across the five datasets
(H1, H2, H3, M1, M2), separately in each immune cell type.

PSI of a junction in a sample = junction reads / reads of the event (LeafCutter cluster), in samples with >= 20 reads
in the event (as in the main analysis). For each junction: one-way ANOVA (scipy.stats.f_oneway) over the five
datasets, the samples of each dataset being one group. Event P = smallest junction P x number of junctions
(Bonferroni within the event, capped at 1); P adjust = Benjamini-Hochberg over the tested events of the cell type.
Tested: every junction of the event has a PSI in at least two samples of each of the five datasets.
unchanged = ANOVA not significant (P adjust > 0.05) and max_range <= 0.1: these events are not tested further.
Also reported: the mean PSI of each dataset for the junction with the smallest P, and the largest range of dataset
mean PSIs over the junctions of the event (max_range).

Output (ensembl116/revision_6ds/): anova_events.tsv (one row per event x cell type), anova_summary.tsv/.xlsx
Usage: anova_events_116.py
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'separation'))
import build_supp_tables_116 as bst  # noqa: E402

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
E = B + 'events/id50_6ds/'
OUT = B + 'revision_6ds/'
DATASETS = ['H1', 'H2', 'H3', 'M1', 'M2']
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
MIN_READS = 20


def main():
    os.makedirs(OUT, exist_ok=True)
    counts = pd.read_csv(E + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(DATASETS)]
    ann = pd.read_csv(B + 'events/id50_sep/classification/events_annotation.tsv', sep='\t',
                      usecols=['cluster', 'symbol', 'event_type'])
    clu_all = pd.Series(counts.index.str.split(':').str[3], index=counts.index)
    rows, summary = [], []
    for cell in CELLS:
        smp = design.index[design['cell'] == cell].intersection(counts.columns)
        tot = counts[smp].groupby(clu_all).transform('sum')
        psi = counts[smp] / tot.where(tot >= MIN_READS)
        groups = {d: design.index[(design['cell'] == cell) & (design['dataset'] == d)].intersection(smp)
                  for d in DATASETS}
        ok = pd.concat([psi[g].notna().sum(axis=1) >= 2 for g in groups.values()], axis=1).all(axis=1)
        ok = ok.groupby(clu_all).transform('all')
        p = psi[ok]
        pv = []
        for j in p.index:
            vals = [p.loc[j, g].dropna().to_numpy() for g in groups.values()]
            if all(np.ptp(v) == 0 for v in vals) and len({float(v[0]) for v in vals}) == 1:
                pv.append(1.0)   # identical PSI everywhere (e.g. all 1): no variance at all
            else:
                pv.append(f_oneway(*vals).pvalue)
        means = pd.DataFrame({d: p[g].mean(axis=1) for d, g in groups.items()})
        j = means.assign(p=np.nan_to_num(pv, nan=1.0), cluster=clu_all[ok].to_numpy(),
                         rng=means.max(axis=1) - means.min(axis=1))
        best = j.loc[j.groupby('cluster')['p'].idxmin()].set_index('cluster')
        ev = pd.DataFrame({'n_junctions': j.groupby('cluster').size(),
                           'p': j.groupby('cluster')['p'].apply(lambda v: min(1.0, v.min() * len(v))),
                           'max_range': j.groupby('cluster')['rng'].max()})
        ev['p_adjust'] = multipletests(ev['p'], method='fdr_bh')[1]
        ev['best_junction'] = best.index.map(lambda c: j.index[(j['cluster'] == c)][j.loc[j['cluster'] == c, 'p'].argmin()])
        for d in DATASETS:
            ev[f'mean_PSI_{d}'] = best[d].round(4)
        ev['anova_significant'] = ev['p_adjust'] <= 0.05
        ev['unchanged'] = ~ev['anova_significant'] & (ev['max_range'] <= 0.1)   # no evidence of change, and small
        ev = ev.reset_index().merge(ann, on='cluster', how='left').assign(cell=cell)
        rows.append(ev)
        n = len(ev)
        ns = ~ev['anova_significant']
        summary.append({'cell': cell, 'tested_events': n, 'significant (P adjust <= 0.05)': int((~ns).sum()),
                        'pct_significant': round(100 * (~ns).mean(), 1), 'not significant': int(ns.sum()),
                        'pct_not_significant': round(100 * ns.mean(), 1),
                        'not significant, max_range <= 0.1': int((ns & (ev['max_range'] <= 0.1)).sum()),
                        'not significant, max_range > 0.1': int((ns & (ev['max_range'] > 0.1)).sum())})
        print(cell, n, 'events', flush=True)
    cols = ['cell', 'cluster', 'symbol', 'event_type', 'n_junctions', 'best_junction'] + \
           [f'mean_PSI_{d}' for d in DATASETS] + ['max_range', 'p', 'p_adjust', 'anova_significant', 'unchanged']
    out = pd.concat(rows)[cols]
    out.to_csv(OUT + 'anova_events.tsv', sep='\t', index=False, float_format='%.6g')
    s = pd.DataFrame(summary)
    s.to_csv(OUT + 'anova_summary.tsv', sep='\t', index=False)
    bst.OUT = OUT
    bst.write('anova_summary', 'One-way ANOVA of junction PSI across the five datasets, per cell type',
              'Per junction: one-way ANOVA of sample PSI (samples with >= 20 reads in the event) with the five datasets '
              '(H1, H2, H3, M1, M2) as groups. Event P = smallest junction P x number of junctions; P adjust = '
              'Benjamini-Hochberg over the tested events of the cell type. Tested: every junction has a PSI in at least '
              'two samples of each dataset. max_range: largest range of dataset mean PSIs over the junctions of the '
              'event. Per-event results: anova_events.tsv.', {'summary': s})
    print(s.to_string(index=False))


if __name__ == '__main__':
    main()
