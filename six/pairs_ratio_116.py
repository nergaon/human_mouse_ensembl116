#!/usr/bin/env python3
"""
Second step after anova_events_116.py: pairwise comparison of the five datasets for the events that are not
'unchanged' (ANOVA significant, or not significant with max_range > 0.1), per immune cell type.

For each event and each of the ten dataset pairs: |dPSI| = largest absolute difference in dataset mean PSI over the
junctions of the event (dataset PSI = mean over samples with >= 20 reads in the event; same samples as the ANOVA).
Six human-mouse pairs (H-M), three human-human pairs (H-H) and one mouse-mouse pair (M-M).
  DS (strict)  : every H-M |dPSI| > 2 x every H-H and M-M |dPSI| (min H-M > 2 x max within)
  DS (mean)    : mean H-M |dPSI| > 2 x mean within-species |dPSI| (H-H and M-M together)
  ratio        : min H-M / max within, and mean H-M / mean within
Also flagged: DS events whose min H-M |dPSI| >= 0.1.

Output (ensembl116/revision_6ds/): pairs_ratio_events.tsv, pairs_ratio_summary.(tsv|xlsx)
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'separation'))
import build_supp_tables_116 as bst  # noqa: E402

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
E = B + 'events/id50_6ds/'
OUT = B + 'revision_6ds/'
DATASETS = ['H1', 'H2', 'H3', 'M1', 'M2']
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
MIN_READS = 20
PAIRS = list(itertools.combinations(DATASETS, 2))
HM = [p for p in PAIRS if p[0][0] != p[1][0]]
WITHIN = [p for p in PAIRS if p[0][0] == p[1][0]]


def main():
    an = pd.read_csv(OUT + 'anova_events.tsv', sep='\t')
    counts = pd.read_csv(E + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(DATASETS)]
    clu_all = pd.Series(counts.index.str.split(':').str[3], index=counts.index)
    rows, summary = [], []
    for cell in CELLS:
        a = an[an['cell'] == cell]
        cand = set(a.loc[~a['unchanged'], 'cluster'])
        smp = design.index[design['cell'] == cell].intersection(counts.columns)
        tot = counts[smp].groupby(clu_all).transform('sum')
        psi = counts[smp] / tot.where(tot >= MIN_READS)
        keep = clu_all.isin(cand)
        means = pd.DataFrame({d: psi.loc[keep, design.index[(design['cell'] == cell) & (design['dataset'] == d)]
                                         .intersection(smp)].mean(axis=1) for d in DATASETS})
        means['cluster'] = clu_all[keep]
        d = pd.DataFrame({f'{x}-{y}': (means[x] - means[y]).abs() for x, y in PAIRS})
        d['cluster'] = means['cluster']
        ev = d.groupby('cluster').max()
        hm = ev[[f'{x}-{y}' for x, y in HM]]
        wi = ev[[f'{x}-{y}' for x, y in WITHIN]]
        ev['min_HM'], ev['max_within'] = hm.min(axis=1), wi.max(axis=1)
        ev['mean_HM'], ev['mean_within'] = hm.mean(axis=1), wi.mean(axis=1)
        ev['ratio_min_HM_max_within'] = ev['min_HM'] / ev['max_within'].replace(0, np.nan)
        ev['ratio_mean'] = ev['mean_HM'] / ev['mean_within'].replace(0, np.nan)
        ev['DS_strict'] = ev['min_HM'] > 2 * ev['max_within']
        ev['DS_mean'] = ev['mean_HM'] > 2 * ev['mean_within']
        ev['DS_strict_and_min_HM_ge_0.1'] = ev['DS_strict'] & (ev['min_HM'] >= 0.1)
        ev = ev.reset_index().assign(cell=cell).merge(a[['cluster', 'symbol', 'event_type']], on='cluster', how='left')
        rows.append(ev)
        summary.append({'cell': cell, 'tested_by_ANOVA': len(a), 'unchanged (not tested further)': int(a['unchanged'].sum()),
                        'tested_pairwise': len(ev), 'DS_strict': int(ev['DS_strict'].sum()),
                        'DS_strict_and_min_HM_ge_0.1': int(ev['DS_strict_and_min_HM_ge_0.1'].sum()),
                        'DS_mean': int(ev['DS_mean'].sum()),
                        'pct_DS_strict_of_tested_by_ANOVA': round(100 * ev['DS_strict'].sum() / len(a), 1),
                        'pct_DS_mean_of_tested_by_ANOVA': round(100 * ev['DS_mean'].sum() / len(a), 1)})
    out = pd.concat(rows)
    lead = ['cell', 'cluster', 'symbol', 'event_type']
    out = out[lead + [c for c in out.columns if c not in lead]]
    out.to_csv(OUT + 'pairs_ratio_events.tsv', sep='\t', index=False, float_format='%.4g')
    s = pd.DataFrame(summary)
    s.to_csv(OUT + 'pairs_ratio_summary.tsv', sep='\t', index=False)
    bst.OUT = OUT
    bst.write('pairs_ratio_summary', 'Pairwise step: human-mouse vs within-species |dPSI| for the events not unchanged',
              '|dPSI| of a pair = largest absolute difference in dataset mean PSI over the junctions of the event. '
              'DS strict: every human-mouse |dPSI| > 2 x every human-human and mouse-mouse |dPSI|. DS mean: mean '
              'human-mouse |dPSI| > 2 x mean within-species |dPSI|. Unchanged (ANOVA not significant and max_range <= '
              '0.1) events are not tested. Per-event values: pairs_ratio_events.tsv.', {'summary': s})
    print(s.to_string(index=False))


if __name__ == '__main__':
    main()
