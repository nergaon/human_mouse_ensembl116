#!/usr/bin/env python3
"""
Calibration of a sample-level Wilcoxon rank-sum test (not part of the manuscript).

For each cell type and each of the ten 3-vs-2 splits of the five datasets, the PSI of each junction in each sample
(samples with >= 20 reads in the event) of the three datasets is compared with that of the two datasets by a
Mann-Whitney U test. Event P = min junction P x number of junctions (Bonferroni within the event), BH across events;
DS = adjusted P <= 0.05 and |difference in mean PSI| >= 0.1 for that junction. Testable events = those of the
five-dataset classification (PSI in all five datasets), so the same events are tested in every split.

Output: ensembl116/events/id50_pooled/summary/wilcoxon_calibration.tsv
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import classify_separation_116 as cs  # noqa: E402

E = cs.BASE + 'events/'
OUT = E + 'id50_pooled/summary/'


def main():
    counts = pd.read_csv(E + 'id50/AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'id50/leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(cs.DATASETS)]
    clu_all = counts.index.str.split(':').str[3]
    rows = []
    for cell in cs.CELLS:
        dpsi, clu = cs.dataset_psi(counts, design, cell)          # testable junctions (PSI in all 5 datasets)
        smp = design.index[design['cell'] == cell].intersection(counts.columns)
        c = counts.loc[dpsi.index, smp]
        tot = counts[smp].groupby(clu_all).transform('sum').loc[dpsi.index]
        psi = c / tot.where(tot >= cs.MIN_READS)
        ds_of = design.loc[smp, 'dataset']
        for a in itertools.combinations(cs.DATASETS, 3):
            ina = ds_of.isin(a).to_numpy()
            x, y = psi.loc[:, ina].to_numpy(), psi.loc[:, ~ina].to_numpy()
            p = np.array([mannwhitneyu(xi[~np.isnan(xi)], yi[~np.isnan(yi)]).pvalue for xi, yi in zip(x, y)])
            d = np.abs(np.nanmean(x, 1) - np.nanmean(y, 1))
            j = pd.DataFrame({'clu': clu, 'p': p, 'big': d >= 0.1})
            ev = j.groupby('clu').agg(p=('p', lambda v: min(1, v.min() * len(v))), n=('p', 'size'))
            # the junction with the smallest P must also have |dPSI| >= 0.1
            best = j.loc[j.groupby('clu')['p'].idxmin()].set_index('clu')['big']
            ev['padj'] = multipletests(ev['p'], method='fdr_bh')[1]
            ds = (ev['padj'] <= 0.05) & best.reindex(ev.index)
            rows.append({'cell': cell, 'split': ' '.join(a) + ' vs ' + ' '.join(d_ for d_ in cs.DATASETS if d_ not in a),
                         'human_mouse': a == cs.HUMAN, 'events': len(ev), 'DS': int(ds.sum()),
                         'pct_DS': round(100 * ds.mean(), 1)})
        print(cell, 'done', flush=True)
    df = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(OUT + 'wilcoxon_calibration.tsv', sep='\t', index=False)
    p = df.pivot(index='split', columns='cell', values='pct_DS')[cs.CELLS]
    hm = df.loc[df['human_mouse'], 'split'].iloc[0]
    print(p.loc[[hm] + [s for s in p.index if s != hm]])
    n = df[~df['human_mouse']].groupby('cell')['DS'].mean()
    h = df[df['human_mouse']].set_index('cell')['DS']
    print('FDR estimate (mean mixed / human-mouse):'); print((n / h).round(2)[cs.CELLS])


if __name__ == '__main__':
    main()
