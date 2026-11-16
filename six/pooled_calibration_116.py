#!/usr/bin/env python3
"""
Summary of the pooled LeafCutter calibration (leafcutter_pooled_design_116.py): for each of the ten 3-vs-2 splits of
the five datasets and each cell type, the number of tested events and of events called differentially spliced
(p.adjust <= 0.05 and max|dPSI| >= 0.1, the single-pair rule), on all tested events and on the events testable in the
five-dataset classification. The human-mouse split is the pooled species test; the nine mixed splits estimate how
many events a pooled test calls from dataset differences alone.

Output: ensembl116/events/id50_pooled/summary/pooled_calibration.tsv
"""
import os

import pandas as pd

E = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/'
LC = E + 'id50_pooled/leafcutter/'
OUT = E + 'id50_pooled/summary/'
CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']
HM = 'H1H2H3_vs_M1M2'


def main():
    st = pd.read_csv(E + 'id50_6ds_anova/classification/event_status_by_cell.tsv', sep='\t',
                     usecols=['cell', 'cluster', 'status_high_low'])
    rows = []
    for split in sorted(os.listdir(LC)):
        if not os.path.isdir(LC + split):
            continue
        for cell in CELLS:
            d = LC + f'{split}/{cell}/'
            sig = pd.read_csv(d + 'leafcutter_ds_cluster_significance.txt', sep='\t')
            sig = sig[sig['status'] == 'Success'].copy()
            sig['cl'] = sig['cluster'].str.split(':').str[1]
            eff = pd.read_csv(d + 'leafcutter_ds_effect_sizes.txt', sep='\t')
            eff['cl'] = eff['intron'].str.split(':').str[3]
            mx = eff.assign(a=eff['deltapsi'].abs()).groupby('cl')['a'].max()
            sig['max_abs_dpsi'] = sig['cl'].map(mx)
            sig['DS'] = (sig['p.adjust'] <= 0.05) & (sig['max_abs_dpsi'] >= 0.1)
            testable = set(st[(st['cell'] == cell) & (st['status_high_low'] != 'not testable')]['cluster'])
            sub = sig[sig['cl'].isin(testable)]
            rows.append({'split': split.replace('_vs_', ' vs '), 'human_mouse': split == HM, 'cell': cell,
                         'tested': len(sig), 'DS': int(sig['DS'].sum()),
                         'pct_DS': round(100 * sig['DS'].mean(), 1),
                         'tested_classification_events': len(sub), 'DS_classification_events': int(sub['DS'].sum()),
                         'pct_DS_classification_events': round(100 * sub['DS'].mean(), 1)})
    df = pd.DataFrame(rows)
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(OUT + 'pooled_calibration.tsv', sep='\t', index=False)
    p = df.pivot(index='split', columns='cell', values='pct_DS')[CELLS]
    p = p.loc[[HM.replace('_vs_', ' vs ')] + [s for s in p.index if s != HM.replace('_vs_', ' vs ')]]
    print('% DS of tested events (p.adjust <= 0.05, max|dPSI| >= 0.1)'); print(p)
    print(df.pivot(index='split', columns='cell', values='DS')[CELLS].loc[p.index])
    n = df[~df['human_mouse']].groupby('cell')['DS'].mean()
    h = df[df['human_mouse']].set_index('cell')['DS']
    print('mean mixed / human-mouse (FDR estimate):'); print((n / h).round(2)[CELLS])
    print(df[['tested', 'tested_classification_events']].describe().loc[['min', 'max']])


if __name__ == '__main__':
    main()
