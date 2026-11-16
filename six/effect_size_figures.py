#!/usr/bin/env python3
"""
Effect-size figures (reviewers 1, 2): human-mouse vs within-species splicing
differences on the SAME events, independent of significance thresholds.

  A  cumulative distribution of event max|dPSI| per cell type
  B  Pearson correlation of group-mean PSI per cell type
  C  group-mean PSI scatter (one cell type) for each comparison

Input: junction_effects.tsv / event_effects.tsv from effect_size_analysis.py
Usage: effect_size_figures.py <effects_dir> [--scatter-cell NK]
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']
CELL_ALIASES = {'BCell': 'NveB', 'Cd8T': 'CD8T'}
# fixed categorical order: slot 1 blue, slot 2 orange, slot 3 aqua
COMPARISONS = {'HH': ('Human-Human', '#2a78d6'), 'HM': ('Human-Mouse', '#eb6834'), 'MM': ('Mouse-Mouse', '#1baf7a')}
INK, MUTED, GRID = '#222222', '#6b6b6b', '#e4e4e4'

plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})


def common_sets(jun, ev):
    """Restrict each cell type to junctions testable in all three comparisons."""
    jun = jun[jun['testable'] & jun['comparison'].isin(COMPARISONS)].copy()
    jun['cell'] = jun['cell'].replace(CELL_ALIASES)
    ev = ev.copy()
    ev['cell'] = ev['cell'].replace(CELL_ALIASES)
    keep_j, keep_e = [], []
    for cell in CELLS:
        t = jun[jun['cell'] == cell]
        common = set.intersection(*[set(t.loc[t['comparison'] == c, 'intron']) for c in COMPARISONS])
        j = t[t['intron'].isin(common)]
        keep_j.append(j)
        clusters = set(j['cluster'])
        keep_e.append(ev[(ev['cell'] == cell) & ev['comparison'].isin(COMPARISONS) & ev['cluster'].isin(clusters)])
    return pd.concat(keep_j), pd.concat(keep_e)


def panel_ecdf(axes, ev):
    for ax, cell in zip(axes, CELLS):
        for comp, (label, color) in COMPARISONS.items():
            x = np.sort(ev.loc[(ev['cell'] == cell) & (ev['comparison'] == comp), 'max_abs_dpsi'].dropna())
            ax.step(x, np.arange(1, len(x) + 1) / len(x), where='post', color=color, lw=1.5, label=label)
        n = ev.loc[(ev['cell'] == cell) & (ev['comparison'] == 'HM')].shape[0]
        ax.axvline(0.1, color=MUTED, lw=0.8, ls=':')
        ax.set_title(f'{cell} (n={n})', fontsize=8, color=INK)
        ax.set_xlim(0, 0.6)
        ax.set_xticks([0, 0.25, 0.5])
        ax.set_ylim(0, 1.01)
        ax.grid(color=GRID, lw=0.5)
        ax.set_xlabel('event max |ΔPSI|')
        if ax is not axes[0]:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel('cumulative fraction\nof events')
    axes[2].legend(frameon=False, ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1.18))


def panel_correlation(ax, jun):
    width = 0.26
    for k, (comp, (label, color)) in enumerate(COMPARISONS.items()):
        r = [stats.pearsonr(*jun.loc[(jun['cell'] == c) & (jun['comparison'] == comp), ['psi1', 'psi2']].T.values)[0]
             for c in CELLS]
        x = np.arange(len(CELLS)) + (k - 1) * width
        ax.bar(x, r, width - 0.04, color=color, label=label)
        for xi, ri in zip(x, r):
            ax.text(xi, ri + 0.01, f'{ri:.2f}', ha='center', va='bottom', fontsize=6, color=INK)
    ax.set_xticks(np.arange(len(CELLS)), CELLS)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Pearson r of junction PSI')
    ax.grid(axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)



def panel_scatter(axes, jun, cell):
    labels = {'HH': ('GSE115736 PSI', 'GSE60424 PSI'), 'HM': ('Human PSI', 'Mouse PSI'),
              'MM': ('GSE116177 PSI', 'GSE180020 PSI')}
    for ax, (comp, (label, color)) in zip(axes, COMPARISONS.items()):
        j = jun[(jun['cell'] == cell) & (jun['comparison'] == comp)]
        ax.hexbin(j['psi1'], j['psi2'], gridsize=30, bins='log', cmap='Greys', mincnt=1, linewidths=0)
        ax.plot([0, 1], [0, 1], color=color, lw=1)
        r = stats.pearsonr(j['psi1'], j['psi2'])[0]
        ax.set_title(f'{label}, {cell}: r = {r:.2f}', fontsize=8, color=INK)
        ax.set_xlabel(labels[comp][0])
        ax.set_ylabel(labels[comp][1])
        ax.set_aspect('equal')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('effects_dir')
    ap.add_argument('--scatter-cell', default='NK')
    args = ap.parse_args()
    jun = pd.read_csv(os.path.join(args.effects_dir, 'junction_effects.tsv'), sep='\t')
    ev = pd.read_csv(os.path.join(args.effects_dir, 'event_effects.tsv'), sep='\t')
    jun, ev = common_sets(jun, ev)

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 5, height_ratios=[1, 0.9, 1.15], hspace=0.65, wspace=0.12)
    panel_ecdf([fig.add_subplot(gs[0, i]) for i in range(5)], ev)
    panel_correlation(fig.add_subplot(gs[1, :]), jun)
    sgs = gs[2, :].subgridspec(1, 3, wspace=0.45)
    panel_scatter([fig.add_subplot(sgs[0, i]) for i in range(3)], jun, args.scatter_cell)
    for letter, y in zip('ABC', (0.9, 0.62, 0.33)):
        fig.text(0.02, y, letter, fontsize=11, fontweight='bold', color=INK)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(args.effects_dir, f'effect_size_figure.{ext}'), dpi=300, bbox_inches='tight')

    # numbers behind the figure, for the text
    rows = []
    for (cell, comp), e in ev.groupby(['cell', 'comparison']):
        j = jun[(jun['cell'] == cell) & (jun['comparison'] == comp)]
        rows.append(dict(cell=cell, comparison=comp, events=len(e), junctions=len(j),
                         median_event_max_abs_dpsi=e['max_abs_dpsi'].median(),
                         frac_events_dpsi_ge_0_1=(e['max_abs_dpsi'] >= 0.1).mean(),
                         pearson=stats.pearsonr(j['psi1'], j['psi2'])[0]))
    pd.DataFrame(rows).to_csv(os.path.join(args.effects_dir, 'effect_size_common_events.tsv'), sep='\t', index=False)
    # one-sided test per cell type: are HM event |dPSI| larger than within-species?
    for cell in CELLS:
        e = ev[ev['cell'] == cell].pivot_table(index='cluster', columns='comparison', values='max_abs_dpsi')
        for other in ('HH', 'MM'):
            p = stats.wilcoxon(e['HM'], e[other], alternative='greater').pvalue
            print(f'{cell}: HM > {other} paired Wilcoxon p = {p:.2e}')


if __name__ == '__main__':
    main()
