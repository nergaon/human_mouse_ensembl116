#!/usr/bin/env python3
"""
Splicing differences of the classified events in every comparison.

For each immune cell type, the events classified as species-divergent or as conserved (high or low confidence, i.e.
from the human-mouse pairs only) are shown with their max|dPSI| in each of the ten comparisons (3 human-human,
1 mouse-mouse, 6 human-mouse) in which they are testable, so that the within-species differences of these events
can be compared with their human-mouse differences.
Figure S3. Output: revision_separation/figures/Figure_S3.pdf/.tif and events/id50_sep/revision/classified_dpsi_by_comparison.tsv
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figure_dpsi_ci_116 import ORDER, boxes
from figures_main_116 import CELL_LABEL, CELLS, ED, KIND_COLORS, KIND_LABEL, kind, label, save

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})
CLASSES = ('divergent', 'conserved')


def main():
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    ev = ev[ev['comparison'].isin(ORDER) & (ev['status_m0.1'] != 'not_testable')]
    st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
    fig, axes = plt.subplots(len(CELLS), 2, figsize=(7.2, 10.5), gridspec_kw={'hspace': 0.75, 'wspace': 0.25})
    rows = []
    for i, cell in enumerate(CELLS):
        for j, cls in enumerate(CLASSES):
            cl = st.loc[(st['cell'] == cell) & (st['status_high_low'] == cls), 'cluster']
            e = ev[(ev['cell'] == cell) & ev['cluster'].isin(cl)]
            data = [e.loc[e['comparison'] == c, 'max_abs_dpsi'].values for c in ORDER]
            for c, d in zip(ORDER, data):
                rows.append(dict(cell=cell, cls=cls, comparison=c, kind=KIND_LABEL[kind(c)], events=len(d),
                                 median_max_abs_dpsi=np.median(d) if len(d) else np.nan,
                                 pct_dpsi_ge_0_1=100 * np.mean(d >= 0.1) if len(d) else np.nan))
            ax = axes[i, j]
            boxes(ax, data, (0, 0.6))
            ax.axhline(0.1, color='#6b6b6b', lw=0.6, ls=':')
            for x, d in enumerate(data):
                ax.text(x, 1.01, f'{np.median(d):.2f}' if len(d) else '', ha='center', va='bottom', fontsize=4.5,
                        transform=ax.get_xaxis_transform())
            ax.set_ylabel('event max |ΔPSI|')
            name = 'Differentially spliced' if cls == 'divergent' else 'Equivalent'
            ax.set_title(f'{CELL_LABEL[cell]}: {name} (n = {len(cl)} events)', fontsize=7, pad=11)
            if i == 0:
                label(ax, 'AB'[j], x=-0.2, y=1.12)
    handles = [plt.Rectangle((0, 0), 1, 1, color=KIND_COLORS[k]) for k in ('HH', 'MM', 'cross')]
    fig.legend(handles, [KIND_LABEL[k] for k in ('HH', 'MM', 'cross')], loc='upper center', ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.93), fontsize=7)
    fig.text(0.02, 0.08, 'Figure S3', fontsize=8)
    save(fig, 'Figure_S3')
    t = pd.DataFrame(rows).round(3)
    t.to_csv(ED + 'revision/classified_dpsi_by_comparison.tsv', sep='\t', index=False)
    print(t.groupby(['cls', 'kind'])['median_max_abs_dpsi'].agg(['min', 'max']).to_string())


if __name__ == '__main__':
    main()
