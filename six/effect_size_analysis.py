#!/usr/bin/env python3
"""
Effect-size and equivalence analysis of orthologous AS events (reviewer requests).

For every comparison (human-mouse, human-human, mouse-mouse, ...) and cell type:
  - per-sample PSI of each junction = junction reads / cluster reads, for samples
    with cluster coverage >= MIN_COVERAGE
  - per-junction group mean PSI, dPSI and a Welch t confidence interval of dPSI
  - per event (LeafCutter cluster): max|dPSI| and the largest |CI bound|
  - status:
      DS            LeafCutter p.adjust <= 0.05 and max|dPSI| >= DS_DPSI
      equivalent_m  every junction's (1 - 2*ALPHA) CI lies within +-m
                    (= two one-sided tests, TOST, at level ALPHA)
      inconclusive  testable but neither of the above
      not_testable  not enough covered samples
  - PSI correlations between the two groups across events

Outputs (OUT_DIR/):
  junction_effects.tsv   one row per comparison x cell type x junction
  event_effects.tsv      one row per comparison x cell type x event, with status
  psi_correlations.tsv   Pearson/Spearman of group-mean PSI per comparison x cell type
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

MIN_COVERAGE = 20      # LeafCutter default min_coverage per sample
ALPHA = 0.05           # TOST level, i.e. 90% two-sided CI
DS_DPSI = 0.1
DS_PADJ = 0.05
MARGINS = (0.025, 0.05, 0.1)

CELL_ALIASES = {'BCell': 'NveB', 'Cd8T': 'CD8T'}


def read_counts(path):
    counts = pd.read_csv(path, sep='\t', index_col=0)
    counts.index.name = 'intron'
    counts['cluster'] = counts.index.str.rsplit(':', n=1).str[1]
    return counts


def min_samples(n1, n2):
    """LeafCutter -g/-i rule used in leafcutter_0.2.9_HN5.sh: 60% of the smaller group, at least 2."""
    return max(2, int(np.floor(0.6 * min(n1, n2))))


def junction_psi(counts, samples):
    sub = counts[samples + ['cluster']]
    totals = sub.groupby('cluster')[samples].transform('sum')
    psi = sub[samples] / totals.where(totals >= MIN_COVERAGE)
    return psi, totals


def welch_ci(x1, x2, level):
    """Mean difference x2 - x1 and its two-sided CI, row-wise; NaN-aware."""
    n1, n2 = x1.notna().sum(axis=1), x2.notna().sum(axis=1)
    m1, m2 = x1.mean(axis=1), x2.mean(axis=1)
    v1, v2 = x1.var(axis=1, ddof=1), x2.var(axis=1, ddof=1)
    se2 = v1 / n1 + v2 / n2
    se = np.sqrt(se2)
    df = se2 ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    # identical PSI in all replicates (se == 0): CI collapses to the point estimate
    df = df.fillna(np.inf).where(se > 0, np.inf)
    t = stats.t.ppf(1 - (1 - level) / 2, df)
    d = m2 - m1
    return d, d - t * se, d + t * se, m1, m2, n1, n2


def analyse(counts, groups, comparison, cell, sig_path):
    g1, g2 = groups.iloc[:, 1].unique()[:2]
    s1 = groups.loc[groups.iloc[:, 1] == g1, groups.columns[0]].tolist()
    s2 = groups.loc[groups.iloc[:, 1] == g2, groups.columns[0]].tolist()
    psi, _ = junction_psi(counts, s1 + s2)
    d, lo, hi, m1, m2, n1, n2 = welch_ci(psi[s1], psi[s2], 1 - 2 * ALPHA)
    need = min_samples(len(s1), len(s2))
    jun = pd.DataFrame({'comparison': comparison, 'cell': cell, 'group1': g1, 'group2': g2,
                        'cluster': counts['cluster'], 'psi1': m1, 'psi2': m2, 'dpsi': d,
                        'ci_low': lo, 'ci_high': hi, 'n1': n1, 'n2': n2})
    jun['testable'] = (jun['n1'] >= need) & (jun['n2'] >= need)

    ev = jun.groupby('cluster').agg(testable=('testable', 'all'), n_junctions=('dpsi', 'size'),
                                    max_abs_dpsi=('dpsi', lambda x: x.abs().max()))
    jun['ci_extreme'] = np.maximum(jun['ci_low'].abs(), jun['ci_high'].abs())
    ev['max_abs_ci'] = jun.groupby('cluster')['ci_extreme'].max()
    # cluster-level PSI vectors for correlations: all junctions of testable clusters

    sig = pd.read_csv(sig_path, sep='\t')
    sig['cluster'] = sig['cluster'].str.split(':').str[1]
    ev = ev.join(sig.set_index('cluster')[['status', 'p.adjust']].rename(columns={'status': 'lc_status'}))
    ev['lc_success'] = ev['lc_status'].eq('Success')
    ds = ev['lc_success'] & (ev['p.adjust'] <= DS_PADJ) & (ev['max_abs_dpsi'] >= DS_DPSI)
    for m in MARGINS:
        st = np.where(ds, 'DS', np.where(ev['testable'] & (ev['max_abs_ci'] < m), 'equivalent', 'inconclusive'))
        ev[f'status_m{m}'] = np.where(ev['testable'] | ds, st, 'not_testable')
    # old classification, for comparison with the manuscript
    old_su = ev['lc_success'] & ~ds & ((ev['p.adjust'] > DS_PADJ) | (ev['max_abs_dpsi'] < 0.05))
    ev['status_old'] = np.where(~ev['lc_success'], 'not_success',
                                np.where(ds, 'DS', np.where(old_su, 'unchanged', 'not_informative')))
    ev.insert(0, 'cell', cell)
    ev.insert(0, 'comparison', comparison)
    return jun.reset_index(), ev.reset_index()


def correlations(jun):
    rows = []
    for (comp, cell), j in jun[jun['testable']].groupby(['comparison', 'cell']):
        x, y = j['psi1'], j['psi2']
        rows.append(dict(comparison=comp, cell=cell, n_junctions=len(j), n_events=j['cluster'].nunique(),
                         pearson=stats.pearsonr(x, y)[0], spearman=stats.spearmanr(x, y)[0],
                         median_abs_dpsi=j['dpsi'].abs().median()))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts', required=True, help='AS_clusters_value table (junction counts)')
    ap.add_argument('--comparison', action='append', required=True,
                    help='name=leafcutter_dir; cell type subdirs hold groups_file.txt and leafcutter_ds_*')
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    counts = read_counts(args.counts)
    juns, evs = [], []
    for spec in args.comparison:
        name, path = spec.split('=', 1)
        for cell in sorted(os.listdir(path)):
            gf = os.path.join(path, cell, 'groups_file.txt')
            sig = os.path.join(path, cell, 'leafcutter_ds_cluster_significance.txt')
            if not (os.path.exists(gf) and os.path.exists(sig)):
                continue
            groups = pd.read_csv(gf, sep='\t', header=None)
            j, e = analyse(counts, groups, name, CELL_ALIASES.get(cell, cell), sig)
            juns.append(j)
            evs.append(e)
            print(name, cell, e[f'status_m0.05'].value_counts().to_dict(), flush=True)
    jun = pd.concat(juns, ignore_index=True)
    ev = pd.concat(evs, ignore_index=True)
    jun.to_csv(os.path.join(args.out_dir, 'junction_effects.tsv'), sep='\t', index=False)
    ev.to_csv(os.path.join(args.out_dir, 'event_effects.tsv'), sep='\t', index=False)
    correlations(jun).to_csv(os.path.join(args.out_dir, 'psi_correlations.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
