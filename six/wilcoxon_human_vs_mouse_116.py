#!/usr/bin/env python3
"""
Per-event Wilcoxon rank-sum (Mann-Whitney) test, all human vs all mouse samples (not part of the manuscript;
separate folder).

For each immune cell type and each event: the PSI of each junction in every sample of H1, H2 and H3 (samples with
>= 20 reads in the event) is compared with its PSI in every sample of M1 and M2 by a two-sided Mann-Whitney U test.
Event P = smallest junction P x number of junctions (Bonferroni within the event, capped at 1); P adjust =
Benjamini-Hochberg over the tested events of the cell type. |max dPSI| = largest absolute difference between the
mean human and the mean mouse sample PSI over the junctions of the event.
Tested events: the events testable in the five-dataset classification (a PSI in every one of the five datasets).

Output: ensembl116/events/id50_pairs_fdr/wilcoxon_human_vs_mouse.(tsv|xlsx), wilcoxon_human_vs_mouse_events.tsv
"""
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
        human = design.loc[smp, 'dataset'].str.startswith('H').to_numpy()
        x, y = psi.loc[:, human].to_numpy(), psi.loc[:, ~human].to_numpy()
        p = np.array([mannwhitneyu(a[~np.isnan(a)], b[~np.isnan(b)], alternative='two-sided').pvalue
                      for a, b in zip(x, y)])
        d = np.abs(np.nanmean(x, axis=1) - np.nanmean(y, axis=1))
        j = pd.DataFrame({'cluster': clu, 'p': p, 'd': d})
        ev = j.groupby('cluster').agg(p=('p', lambda v: min(1.0, v.min() * len(v))), n_junctions=('p', 'size'),
                                      max_abs_dpsi=('d', 'max'))
        ev['p_adjust'] = multipletests(ev['p'], method='fdr_bh')[1]
        sig = ev['p_adjust'] <= 0.05
        sd = sig & (ev['max_abs_dpsi'] >= 0.1)
        rows.append({'Cell type': CELL_NAME[cell], 'Human samples': int(human.sum()),
                     'Mouse samples': int((~human).sum()), 'Tested events': len(ev),
                     'Significant (P adjust <= 0.05)': int(sig.sum()), '% significant': round(100 * sig.mean(), 1),
                     'Significant and |max dPSI| >= 0.1': int(sd.sum()),
                     '% significant and |max dPSI| >= 0.1': round(100 * sd.mean(), 1)})
        events.append(ev.reset_index().assign(cell=cell))
    t = pd.DataFrame(rows)
    t.to_csv(OUT + 'wilcoxon_human_vs_mouse.tsv', sep='\t', index=False)
    pd.concat(events).to_csv(OUT + 'wilcoxon_human_vs_mouse_events.tsv', sep='\t', index=False, float_format='%.6g')
    bst.OUT = OUT
    bst.write('wilcoxon_human_vs_mouse', 'Per-event Wilcoxon rank-sum test: all human vs all mouse samples',
              'Per junction: two-sided Mann-Whitney test of sample PSI (samples with >= 20 reads in the event), all '
              'samples of H1, H2 and H3 vs all samples of M1 and M2 of the cell type. Event P = smallest junction P x '
              'number of junctions; P adjust = Benjamini-Hochberg over the tested events. Tested: the events testable '
              'in the five-dataset classification. |max dPSI|: largest difference in mean sample PSI, human vs mouse.',
              {'by_cell_type': t})
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()
