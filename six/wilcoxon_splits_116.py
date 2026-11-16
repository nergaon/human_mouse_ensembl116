#!/usr/bin/env python3
"""
Per-event Wilcoxon rank-sum (Mann-Whitney) test for each of the ten 3-vs-2 splits of the five datasets (not part of
the manuscript; separate folder). The human-mouse split (H1 H2 H3 vs M1 M2) is the species test; the nine mixed splits
show how many events the same test calls when the groups mix the species (FDR estimate = mean of the nine / the
human-mouse split).

For each immune cell type, split and event: the PSI of each junction in every sample of the three datasets (samples
with >= 20 reads in the event) is compared with its PSI in every sample of the two datasets by a two-sided
Mann-Whitney U test.
Event P = smallest junction P x number of junctions (Bonferroni within the event, capped at 1); P adjust =
Benjamini-Hochberg over the tested events of the cell type. |max dPSI| = largest absolute difference between the
mean sample PSI of the two groups over the junctions of the event.
Tested events: the events testable in the five-dataset classification (a PSI in every one of the five datasets).

Output: ensembl116/events/id50_pairs_fdr/wilcoxon_splits.(tsv|xlsx), wilcoxon_splits_events.tsv
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
import build_supp_tables_116 as bst  # noqa: E402

E = cs.BASE + 'events/id50_6ds/'
OUT = cs.BASE + 'events/id50_pairs_fdr/'
CELL_NAME = {'CD4T': 'CD4 T', 'CD8T': 'CD8 T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Monocytes'}


def main():
    os.makedirs(OUT, exist_ok=True)
    counts = pd.read_csv(E + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(cs.DATASETS)]
    clu_all = counts.index.str.split(':').str[3]
    rows, events = [], []
    for cell in cs.CELLS:
        dpsi, clu = cs.dataset_psi(counts, design, cell)       # testable junctions (PSI in all 5 datasets)
        smp = design.index[design['cell'] == cell].intersection(counts.columns)
        tot = counts[smp].groupby(clu_all).transform('sum').loc[dpsi.index]
        psi = counts.loc[dpsi.index, smp] / tot.where(tot >= cs.MIN_READS)
        dsl = design.loc[smp, 'dataset']
        for grp in itertools.combinations(cs.DATASETS, 3):
            other = tuple(d_ for d_ in cs.DATASETS if d_ not in grp)
            name = ' '.join(grp) + ' vs ' + ' '.join(other)
            ina = dsl.isin(grp).to_numpy()
            x, y = psi.loc[:, ina].to_numpy(), psi.loc[:, ~ina].to_numpy()
            p = np.array([mannwhitneyu(a_[~np.isnan(a_)], b_[~np.isnan(b_)], alternative='two-sided').pvalue
                          for a_, b_ in zip(x, y)])
            d = np.abs(np.nanmean(x, axis=1) - np.nanmean(y, axis=1))
            j = pd.DataFrame({'cluster': clu, 'p': p, 'd': d})
            ev = j.groupby('cluster').agg(p=('p', lambda v: min(1.0, v.min() * len(v))), n_junctions=('p', 'size'),
                                          max_abs_dpsi=('d', 'max'))
            ev['p_adjust'] = multipletests(ev['p'], method='fdr_bh')[1]
            sig = ev['p_adjust'] <= 0.05
            sd = sig & (ev['max_abs_dpsi'] >= 0.1)
            rows.append({'Cell type': CELL_NAME[cell], 'Split': name, 'Human vs mouse': grp == cs.HUMAN,
                         'Samples, group of 3': int(ina.sum()), 'Samples, group of 2': int((~ina).sum()),
                         'Tested events': len(ev),
                         'Significant (P adjust <= 0.05)': int(sig.sum()), '% significant': round(100 * sig.mean(), 1),
                         'Significant and |max dPSI| >= 0.1': int(sd.sum()),
                         '% significant and |max dPSI| >= 0.1': round(100 * sd.mean(), 1)})
            events.append(ev.reset_index().assign(cell=cell, split=name))
        print(cell, 'done', flush=True)
    t = pd.DataFrame(rows)
    fdr = []
    for cell, g in t.groupby('Cell type', sort=False):
        hm, mix = g[g['Human vs mouse']], g[~g['Human vs mouse']]
        r = {'Cell type': cell}
        for k in ('Significant (P adjust <= 0.05)', 'Significant and |max dPSI| >= 0.1'):
            r[f'{k}: human vs mouse'] = int(hm[k].iloc[0])
            r[f'{k}: mean of 9 mixed'] = round(mix[k].mean(), 1)
            r[f'{k}: max of 9 mixed'] = int(mix[k].max())
            r[f'{k}: FDR'] = round(mix[k].mean() / hm[k].iloc[0], 3)
        fdr.append(r)
    fdr = pd.DataFrame(fdr)
    t.to_csv(OUT + 'wilcoxon_splits.tsv', sep='\t', index=False)
    fdr.to_csv(OUT + 'wilcoxon_splits_fdr.tsv', sep='\t', index=False)
    pd.concat(events).to_csv(OUT + 'wilcoxon_splits_events.tsv', sep='\t', index=False, float_format='%.6g')
    bst.OUT = OUT
    bst.write('wilcoxon_splits', 'Per-event Wilcoxon rank-sum test for the ten 3-vs-2 splits of the five datasets',
              'Per junction: two-sided Mann-Whitney test of sample PSI (samples with >= 20 reads in the event), all '
              'samples of the three datasets vs all samples of the two datasets of the cell type. Event P = smallest '
              'junction P x number of junctions; P adjust = Benjamini-Hochberg over the tested events. Tested: the '
              'events testable in the five-dataset classification (the same events in every split). FDR = mean of the '
              'nine mixed splits / human vs mouse split.', {'splits': t, 'FDR': fdr})
    pd.set_option('display.width', 250)
    print(t.to_string(index=False))
    print(fdr.to_string(index=False))


if __name__ == '__main__':
    main()
