#!/usr/bin/env python3
"""
Candidate figure: number of events separated by each of the ten 3-vs-2 splits of the five datasets.

For each cell type and each split of H1, H2, H3, M1, M2 into a group of three and a group of two, an event is
separated when, for at least one junction, the PSI of every dataset of one group differs from the PSI of every dataset
of the other group by at least a gap g (dataset PSI and testable events exactly as in classify_separation_116.py).
The human-mouse split (H1 H2 H3 | M1 M2) gives the differentially spliced events of Table 1; the nine other splits
mix the species and estimate the events separated by dataset differences alone (their mean = expected false
positives, FDR = mean / human-mouse count). 

Output: ensembl116/revision_separation/candidate_figures/Figure_split_counts.(pdf|png), split_counts.tsv
"""
import itertools
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import classify_separation_116 as cs  # noqa: E402

matplotlib.rcParams['font.family'] = 'Liberation Sans'
matplotlib.rcParams['pdf.fonttype'] = 42

ED = cs.BASE + 'events/id50/'
OUT = cs.BASE + 'revision_separation/candidate_figures/'
CELL_LABEL = {'CD4T': 'T4', 'CD8T': 'T8', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo'}


CELL_COLOR = {'CD4T': '#1b9e77', 'CD8T': '#7570b3', 'NveB': '#e7298a', 'NK': '#66a61e', 'Mono': '#e6ab02'}


def split_counts():
    counts = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(ED + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    splits = [cs.HUMAN] + [s for s in itertools.combinations(cs.DATASETS, 3) if s != cs.HUMAN]
    rows = []
    for cell in cs.CELLS:
        psi, clu = cs.dataset_psi(counts, design, cell)
        for s in splits:
            g = cs.gap(psi, s).groupby(clu).max()
            other = tuple(d for d in cs.DATASETS if d not in s)
            row = {'cell': cell, 'split': ' '.join(s) + ' | ' + ' '.join(other), 'human_mouse': s == cs.HUMAN,
                   'testable_events': len(g)}
            for lab, cut in (('0.1', cs.GAP_LOW), ('0.2', cs.GAP_HIGH)):
                sep = g.index[g >= cut]
                row[f'events_gap_ge_{lab}'] = len(sep)
            rows.append(row)
    return pd.DataFrame(rows)


def plot(df):
    splits = df['split'].drop_duplicates().tolist()
    x = np.arange(len(splits))
    w = 0.16
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.2), sharex=True)
    for ax, lab, title in ((axes[0], '0.1', 'Gap ≥ 0.1 (high + low confidence)'),
                           (axes[1], '0.2', 'Gap ≥ 0.2 (high confidence)')):
        ax.axvspan(-0.5, 0.5, color='#eef3f9', zorder=0)
        for k, cell in enumerate(cs.CELLS):
            d = df[df['cell'] == cell].set_index('split').loc[splits]
            v = d[f'events_gap_ge_{lab}'].to_numpy()
            ax.bar(x + (k - 2) * w, v, width=w, color=CELL_COLOR[cell], edgecolor='white', lw=0.3, zorder=2,
                   label=f'{CELL_LABEL[cell]} ({d["testable_events"].iloc[0]:,} testable events)')
            ax.text(x[0] + (k - 2) * w, v[0], str(v[0]), ha='center', va='bottom', fontsize=6, rotation=90)
        null = df[~df['human_mouse']].groupby('cell')[f'events_gap_ge_{lab}'].mean()
        ax.axhline(null.mean(), color='black', ls='--', lw=0.8, zorder=3)
        ax.text(0.55, null.mean(), f' mean of the 9 mixed splits: {null.mean():.0f}', ha='left',
                va='bottom', fontsize=7, zorder=4, bbox=dict(facecolor='white', edgecolor='none', pad=0.5))
        ax.set_ylabel('Separated events', fontsize=9)
        ax.set_title(title, fontsize=10, loc='left')
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylim(0, df[f'events_gap_ge_{lab}'].max() * 1.15)
    axes[0].legend(fontsize=7.5, frameon=False, ncol=3, loc='upper right')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([s.replace(' | ', '\nvs\n') for s in splits], fontsize=8, linespacing=1.0)
    axes[1].set_xlim(-0.5, len(splits) - 0.5)
    axes[1].set_xlabel('Split of the five datasets into a group of three and a group of two '
                       '(shaded: human vs mouse; the other nine mix the species)', fontsize=8.5)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(OUT + f'Figure_split_counts.{ext}', dpi=300)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = split_counts()
    df.to_csv(OUT + 'split_counts.tsv', sep='\t', index=False)
    plot(df)
    for v in ('events_gap_ge_0.1', 'events_gap_ge_0.2'):
        print(v); print(df.pivot(index='split', columns='cell', values=v).loc[df['split'].drop_duplicates()])


if __name__ == '__main__':
    main()
