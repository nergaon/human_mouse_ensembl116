#!/usr/bin/env python3
"""
Distribution of the splicing differences of every comparison, per cell type (reviewers: show that dPSI is larger
between human and mouse than within species).

Events: per cell type, the events testable in all ten comparisons (3 human-human, 1 mouse-mouse, 6 human-mouse),
so every comparison is measured on the same events.
  left   event max|dPSI| (largest |dPSI| of the junctions of the event); dotted line: DS threshold 0.1
  right  largest |bound| of the 90% CI of dPSI over the junctions of the event; an event is equivalent at margin m
         when this value is < m. Lines at 0.05 and 0.1; numbers above each box: % of events < 0.05 / < 0.1
Figure S1. Output: revision_files/figures/Figure_S1.pdf/.tif and events/id50/revision/dpsi_ci_by_comparison.tsv
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures_main_116 import CELL_LABEL, CELLS, ED, KIND_COLORS, KIND_LABEL, OUT, kind, label, save

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})
ORDER = ['H1_H2', 'H1_H3', 'H2_H3', 'M1_M2', 'H1_M1', 'H1_M2', 'H2_M1', 'H2_M2', 'H3_M1', 'H3_M2']


def boxes(ax, data, ylim):
    bp = ax.boxplot(data, positions=range(len(ORDER)), widths=0.6, whis=(5, 95), showfliers=False,
                    patch_artist=True, medianprops=dict(color='black', lw=0.8),
                    whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
    for patch, comp in zip(bp['boxes'], ORDER):
        patch.set_facecolor(KIND_COLORS[kind(comp)])
        patch.set_alpha(0.85)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([c.replace('_', '-') for c in ORDER], rotation=90, fontsize=6)
    ax.set_xlim(-0.6, len(ORDER) - 0.4)
    ax.set_ylim(*ylim)


def main():
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    ev = ev[ev['comparison'].isin(ORDER)]
    fig, axes = plt.subplots(len(CELLS), 2, figsize=(7.2, 10.5), gridspec_kw={'hspace': 0.75, 'wspace': 0.25})
    rows = []
    for i, cell in enumerate(CELLS):
        ec = ev[(ev['cell'] == cell) & ev['testable']]
        n_comp = ec.groupby('cluster')['comparison'].nunique()
        common = n_comp[n_comp == len(ORDER)].index
        ec = ec[ec['cluster'].isin(common)].set_index('comparison')
        dpsi = [ec.loc[c, 'max_abs_dpsi'].values for c in ORDER]
        ci = [ec.loc[c, 'max_abs_ci'].values for c in ORDER]
        for c, d, b in zip(ORDER, dpsi, ci):
            rows.append(dict(cell=cell, comparison=c, kind=KIND_LABEL[kind(c)], events=len(d),
                             median_max_abs_dpsi=np.median(d), pct_dpsi_ge_0_1=100 * np.mean(d >= 0.1),
                             median_ci_bound=np.median(b), pct_ci_lt_0_05=100 * np.mean(b < 0.05),
                             pct_ci_lt_0_1=100 * np.mean(b < 0.1)))

        ax = axes[i, 0]
        boxes(ax, dpsi, (0, 0.6))
        ax.axhline(0.1, color='#6b6b6b', lw=0.6, ls=':')
        for x, d in enumerate(dpsi):
            ax.text(x, 1.01, f'{100 * np.mean(d >= 0.1):.0f}', ha='center', va='bottom', fontsize=5,
                    transform=ax.get_xaxis_transform())
        ax.set_ylabel('event max |ΔPSI|')
        ax.set_title(f'{CELL_LABEL[cell]}: ΔPSI (n = {len(common)} events)', fontsize=7, pad=11)

        ax = axes[i, 1]
        boxes(ax, ci, (0, 0.9))
        for m, ls in ((0.05, '--'), (0.1, ':')):
            ax.axhline(m, color='#6b6b6b', lw=0.6, ls=ls)
        for x, b in enumerate(ci):
            ax.text(x, 1.01, f'{100 * np.mean(b < 0.05):.0f}/{100 * np.mean(b < 0.1):.0f}', ha='center',
                    va='bottom', fontsize=4.5, transform=ax.get_xaxis_transform())
        ax.set_ylabel('largest |90% CI bound| of ΔPSI')
        ax.set_title(f'{CELL_LABEL[cell]}: 90% CI of ΔPSI', fontsize=7, pad=11)
        if i == 0:
            label(axes[i, 0], 'A', x=-0.2, y=1.12)
            label(axes[i, 1], 'B', x=-0.2, y=1.12)
    handles = [plt.Rectangle((0, 0), 1, 1, color=KIND_COLORS[k]) for k in ('HH', 'MM', 'cross')]
    fig.legend(handles, [KIND_LABEL[k] for k in ('HH', 'MM', 'cross')], loc='upper center', ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.93), fontsize=7)
    fig.text(0.02, 0.08, 'Figure S1', fontsize=8)
    save(fig, 'Figure_S1')
    t = pd.DataFrame(rows).round(3)
    t.to_csv(ED + 'revision/dpsi_ci_by_comparison.tsv', sep='\t', index=False)
    print(t.groupby('kind')[['median_max_abs_dpsi', 'pct_dpsi_ge_0_1', 'median_ci_bound', 'pct_ci_lt_0_05',
                             'pct_ci_lt_0_1']].agg(['min', 'max']).to_string())


if __name__ == '__main__':
    main()
