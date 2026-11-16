#!/usr/bin/env python3
"""
Per-event Wilcoxon rank-sum (Mann-Whitney) test for each single human-mouse dataset pair (not part of the
manuscript; separate folder).

For each immune cell type and each of the six human-mouse pairs (H1, H2 or H3 vs M1 or M2): the PSI of each junction
in each sample (samples with >= 20 reads in the event) of the human dataset is compared with that of the mouse dataset
by a two-sided Mann-Whitney U test. Event P = smallest junction P x number of junctions (Bonferroni within the
event, capped at 1); adjusted P = Benjamini-Hochberg over the events tested in the pair. An event is tested when
every junction has a PSI in at least two samples of each dataset. |max dPSI| = largest absolute difference in mean
sample PSI between the two datasets over the junctions of the event.
min_possible_P: the smallest two-sided P the test can give with the pair's sample numbers (exact, no ties).

Output: ensembl116/events/id50_pairs_fdr/wilcoxon_pairs.(tsv|xlsx), wilcoxon_pairs_events.tsv (per event)
"""
import os
import sys
from math import comb

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import classify_separation_116 as cs  # noqa: E402
import build_supp_tables_116 as bst  # noqa: E402

E = cs.BASE + 'events/id50_6ds/'
OUT = cs.BASE + 'events/id50_pairs_fdr/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
HM = [(h, m) for h in ('H1', 'H2', 'H3') for m in ('M1', 'M2')]


def main():
    os.makedirs(OUT, exist_ok=True)
    counts = pd.read_csv(E + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    clu_all = pd.Series(counts.index.str.split(':').str[3], index=counts.index)
    rows, events = [], []
    for cell in CELLS:
        smp = design.index[(design['cell'] == cell)].intersection(counts.columns)
        tot = counts[smp].groupby(clu_all).transform('sum')
        psi = counts[smp] / tot.where(tot >= cs.MIN_READS)
        for h, m in HM:
            sh = design.index[(design['cell'] == cell) & (design['dataset'] == h)].intersection(smp)
            sm = design.index[(design['cell'] == cell) & (design['dataset'] == m)].intersection(smp)
            x, y = psi[sh], psi[sm]
            ok = (x.notna().sum(axis=1) >= 2) & (y.notna().sum(axis=1) >= 2)
            ok = ok.groupby(clu_all).transform('all')
            xs, ys, clu = x[ok].to_numpy(), y[ok].to_numpy(), clu_all[ok].to_numpy()
            p = np.array([mannwhitneyu(a[~np.isnan(a)], b[~np.isnan(b)], alternative='two-sided').pvalue
                          for a, b in zip(xs, ys)])
            d = np.abs(np.nanmean(xs, axis=1) - np.nanmean(ys, axis=1))
            j = pd.DataFrame({'cluster': clu, 'p': p, 'd': d})
            ev = j.groupby('cluster').agg(p=('p', lambda v: min(1.0, v.min() * len(v))), n_junctions=('p', 'size'),
                                          max_abs_dpsi=('d', 'max'))
            ev['p_adjust'] = multipletests(ev['p'], method='fdr_bh')[1]
            sig = ev['p_adjust'] <= 0.05
            sd = sig & (ev['max_abs_dpsi'] >= 0.1)
            n1, n2 = len(sh), len(sm)
            rows.append({'cell': cell, 'pair': f'{h}-{m}', 'human_samples': n1, 'mouse_samples': n2,
                         'min_possible_P': round(2 / comb(n1 + n2, n1), 4), 'tested': len(ev),
                         'smallest_P': ev['p'].min(), 'smallest_P_adjust': ev['p_adjust'].min(),
                         'sig_P_adjust_le_0.05': int(sig.sum()), 'pct_sig': round(100 * sig.mean(), 1),
                         'sig_and_abs_max_dPSI_ge_0.1': int(sd.sum()), 'pct_sig_and_dPSI': round(100 * sd.mean(), 1),
                         'nominal_P_le_0.05': int((ev['p'] <= 0.05).sum())})
            events.append(ev.reset_index().assign(cell=cell, pair=f'{h}-{m}'))
        print(cell, 'done', flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(OUT + 'wilcoxon_pairs.tsv', sep='\t', index=False)
    pd.concat(events).to_csv(OUT + 'wilcoxon_pairs_events.tsv', sep='\t', index=False, float_format='%.6g')
    bst.OUT = OUT
    bst.write('wilcoxon_pairs', 'Single human-mouse dataset pairs: per-event Wilcoxon rank-sum test',
              'Per junction: two-sided Mann-Whitney test of sample PSI (samples with >= 20 reads in the event), human '
              'dataset vs mouse dataset. Event P = smallest junction P x number of junctions; P adjust = '
              'Benjamini-Hochberg over the events tested in the pair. Significant: P adjust <= 0.05; and also '
              '|max dPSI| >= 0.1 (difference in mean sample PSI). min_possible_P: smallest two-sided P possible with '
              'the sample numbers of the pair (before the correction for junctions and events).', {'pairs': t})
    pd.set_option('display.width', 250)
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()
