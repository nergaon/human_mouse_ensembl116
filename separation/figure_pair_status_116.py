#!/usr/bin/env python3
"""
Candidate figure: status of the testable events in each of the ten dataset pairs separately (Figure 2A-C shows the
median over the pairs of each kind). Single-pair test (Figure S1C): LeafCutter p.adjust <= 0.05 and max|dPSI| >= 0.1
-> differentially spliced; 90% CI of dPSI of every junction within +-0.1 -> equivalent; otherwise inconclusive.
Testable = enough samples with >= 20 reads in both datasets (per pair, so the testable events differ between pairs).

Output: ensembl116/revision_separation/candidate_figures/Figure_pair_status.(pdf|png), pair_status.tsv
"""
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})
ED = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50_sep/'
OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_separation/candidate_figures/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_LABEL = {'CD4T': 'T4', 'CD8T': 'T8', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo'}
STATUS = 'status_m0.1'
COLORS = {'DS': '#74c476', 'equivalent': '#9ecae1', 'inconclusive': '#d9d9d9'}   # as in Figure 2
LABEL = {'DS': 'Differentially spliced', 'equivalent': 'Equivalent', 'inconclusive': 'Inconclusive'}
PAIRS = ['H1_H2', 'H1_H3', 'H2_H3', 'M1_M2', 'H1_M1', 'H1_M2', 'H2_M1', 'H2_M2', 'H3_M1', 'H3_M2']
KIND = {'H': 'Human-Human', 'M': 'Mouse-Mouse', 'X': 'Human-Mouse'}
KIND_COLOR = {'H': '#2a78d6', 'M': '#1baf7a', 'X': '#eb6834'}   # KIND_COLORS of figures_main_116.py


def kind(p):
    a, b = p.split('_')
    return 'X' if a[0] != b[0] else a[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t', usecols=['comparison', 'cell', 'cluster', STATUS])
    t = ev[ev['comparison'].isin(PAIRS) & (ev[STATUS] != 'not_testable')]
    n = pd.crosstab([t['comparison'], t['cell']], t[STATUS]).reindex(columns=list(COLORS), fill_value=0)
    pct = n.div(n.sum(axis=1), axis=0) * 100
    out = n.add_prefix('n_').join(pct.round(1).add_prefix('pct_'))
    out['testable'] = n.sum(axis=1)
    out.to_csv(OUT + 'pair_status.tsv', sep='\t')

    fig, axes = plt.subplots(2, 5, figsize=(7.2, 5.0), sharey=True)
    for ax, p in zip(axes.flat, PAIRS):
        d, c = pct.loc[p].reindex(CELLS), n.loc[p].reindex(CELLS)
        x = [CELL_LABEL[k] for k in CELLS]
        bottom = np.zeros(len(CELLS))
        for s, col in COLORS.items():
            ax.bar(x, d[s], bottom=bottom, color=col, width=0.75, edgecolor='white', linewidth=0.5)
            for i, (b, v, k) in enumerate(zip(bottom, d[s], c[s])):
                if v >= 7:
                    ax.text(i, b + v / 2, f'{k}', ha='center', va='center', fontsize=4.8, color='#333333')
            bottom += d[s].to_numpy()
        a, b = p.split('_')
        ax.set_title(f'{a} vs {b}', fontsize=8, color=KIND_COLOR[kind(p)], fontweight='bold')
        ax.set_ylim(0, 100)
        ax.tick_params(axis='x', length=0)
    for ax in axes[:, 0]:
        ax.set_ylabel('% of testable events')
    handles = [Patch(color=c, label=LABEL[s]) for s, c in COLORS.items()]
    handles += [Patch(color=KIND_COLOR[k], label=f'{KIND[k]} pair (title)') for k in 'HMX']
    fig.legend(handles=handles, loc='lower center', ncol=6, fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0),
               handlelength=1.2, columnspacing=1.0)
    fig.tight_layout(rect=(0, 0.05, 1, 1), h_pad=1.5)
    for ext in ('pdf', 'png'):
        fig.savefig(OUT + f'Figure_pair_status.{ext}', dpi=300)
    plt.close(fig)
    print(pct.round(1).unstack('cell').swaplevel(axis=1)['CD4T'])
    print(pct['DS'].groupby(lambda i: kind(i[0])).agg(['min', 'max']).round(1))
    print(pct['equivalent'].groupby(lambda i: kind(i[0])).agg(['min', 'max']).round(1))


if __name__ == '__main__':
    main()
