#!/usr/bin/env python3
"""
Candidate figure: LeafCutter significance of the same events in the human-mouse pair H2-M1 (x) and in a
within-species pair (y), per cell type.

Row 1: x = -log10 p.adjust H2 vs M1, y = -log10 p.adjust H1 vs H2.
Row 2: x = -log10 p.adjust H2 vs M1, y = -log10 p.adjust M1 vs M2.
Events = LeafCutter clusters tested successfully in both pairs of the panel. Dashed lines: p.adjust = 0.05;
rho = Spearman correlation of the two -log10 p.adjust values over the events of the panel. Both axes end at one
common limit, the smallest value that contains 90% of the points of all panels (both coordinates); points beyond it
are drawn on the edge (triangles).

Row 3: estimated FDR of H2 vs M1 at each -log10 p.adjust cutoff = number of events passing the cutoff in the
within-species pair / number passing in H2 vs M1 (same events); the within-species pair estimates how many events
pass from dataset differences alone.

Output: ensembl116/revision_6ds/candidate_figures/Figure_padj_H2M1_vs_within.(pdf|png), padj_pairs.tsv,
        padj_pairs_fdr.tsv
"""
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})
B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
LC = B + 'events/id50_6ds/leafcutter/original/'
OUT = B + 'revision_6ds/candidate_figures/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_NAME = {'CD4T': 'CD4 T cells', 'CD8T': 'CD8 T cells', 'NveB': 'B cells', 'NK': 'NK cells', 'Mono': 'Monocytes'}
X = 'H2_M1'
ROWS = [('H1_H2', '#2a78d6'), ('M1_M2', '#1baf7a')]   # KIND_COLORS of figures_main_116.py
CUT = -np.log10(0.05)


def padj(pair, cell):
    d = pd.read_csv(LC + f'{pair}/{cell}/leafcutter_ds_cluster_significance.txt', sep='\t')
    d = d[d['status'] == 'Success']
    return pd.Series(-np.log10(d['p.adjust'].clip(lower=1e-300)).to_numpy(), index=d['cluster'].str.split(':').str[1])


def main():
    os.makedirs(OUT, exist_ok=True)
    data = {(ypair, cell): pd.concat([padj(X, cell).rename('x'), padj(ypair, cell).rename('y')], axis=1,
                                     join='inner') for ypair, _ in ROWS for cell in CELLS}
    allp = pd.concat(data.values())
    CAP = float(np.ceil(np.quantile(np.maximum(allp['x'], allp['y']), 0.9)))
    print('axis limit (90% of points):', CAP)
    fig, axes = plt.subplots(3, 5, figsize=(10, 6.9))
    rows = []
    for r, (ypair, col) in enumerate(ROWS):
        for c, cell in enumerate(CELLS):
            ax = axes[r, c]
            d = data[(ypair, cell)]
            x, y = d['x'].clip(upper=CAP), d['y'].clip(upper=CAP)
            capped = (d['x'] > CAP) | (d['y'] > CAP)
            ax.scatter(x[~capped], y[~capped], s=3, color=col, alpha=0.35, lw=0, rasterized=True)
            ax.scatter(x[capped], y[capped], s=6, color=col, alpha=0.6, lw=0, marker='^', rasterized=True)
            ax.axvline(CUT, color='#555555', ls='--', lw=0.6)
            ax.axhline(CUT, color='#555555', ls='--', lw=0.6)
            ax.plot([0, CAP], [0, CAP], color='#bbbbbb', lw=0.6, zorder=0)
            q = {'x only': int(((d['x'] >= CUT) & (d['y'] < CUT)).sum()),
                 'y only': int(((d['x'] < CUT) & (d['y'] >= CUT)).sum()),
                 'both': int(((d['x'] >= CUT) & (d['y'] >= CUT)).sum()),
                 'neither': int(((d['x'] < CUT) & (d['y'] < CUT)).sum())}
            rho = spearmanr(d['x'], d['y']).statistic
            ax.set_title(f'{CELL_NAME[cell]}\n{len(d):,} events, rho = {rho:.2f}', fontsize=7)
            ax.set_xlim(-0.3, CAP * 1.03)
            ax.set_ylim(-0.3, CAP * 1.03)
            if c == 0:
                a, b = ypair.split('_')
                ax.set_ylabel(f'-log10 p.adjust, {a} vs {b}')
            ax.set_xlabel('-log10 p.adjust, H2 vs M1', fontsize=6.5)
            rows.append({'y_pair': ypair, 'cell': cell, 'events': len(d), 'spearman': round(rho, 3),
                         'sig_H2_M1_only': q['x only'], f'sig_y_only': q['y only'], 'sig_both': q['both'],
                         'sig_neither': q['neither'],
                         'pct_sig_H2_M1': round(100 * (d['x'] >= CUT).mean(), 1),
                         'pct_sig_y': round(100 * (d['y'] >= CUT).mean(), 1)})
    # row 3: estimated FDR at each cutoff = events passing in the within-species pair / events passing in H2 vs M1
    grid = np.round(np.arange(0, CAP + 0.01, 0.1), 2)
    fdr = []
    for c, cell in enumerate(CELLS):
        ax = axes[2, c]
        for ypair, col in ROWS:
            d = data[(ypair, cell)]
            n_x = np.array([(d['x'] >= t).sum() for t in grid])
            n_y = np.array([(d['y'] >= t).sum() for t in grid])
            ratio = np.where(n_x > 0, n_y / np.maximum(n_x, 1), np.nan)
            a, b = ypair.split('_')
            ax.plot(grid, ratio, color=col, lw=1.2, label=f'{a} vs {b} / H2 vs M1')
            fdr += [{'cell': cell, 'within_pair': ypair, 'neg_log10_padj_cutoff': t, 'events_H2_M1': int(nx),
                     'events_within': int(ny), 'estimated_FDR': round(rt, 4)} for t, nx, ny, rt in
                    zip(grid, n_x, n_y, ratio)]
        ax.axvline(CUT, color='#555555', ls='--', lw=0.6)
        ax.set_xlim(0, CAP)
        ax.set_ylim(0, 1)
        ax.set_xlabel('-log10 p.adjust cutoff', fontsize=6.5)
        if c == 0:
            ax.set_ylabel('Estimated FDR\n(within-species / H2 vs M1)')
            ax.legend(fontsize=5.5, frameon=False, loc='upper right')
    fdr = pd.DataFrame(fdr)
    fdr.to_csv(OUT + 'padj_pairs_fdr.tsv', sep='\t', index=False)
    at = fdr[fdr['neg_log10_padj_cutoff'].isin([1.3, 2.0, 5.0, 10.0])]
    print(at.pivot_table(index=['within_pair', 'neg_log10_padj_cutoff'], columns='cell',
                         values='estimated_FDR')[CELLS].round(2))
    fig.text(0.5, 0.005, f'Rows 1-2: dashed lines, p.adjust = 0.05; rho, Spearman correlation; axes end at {CAP:.0f}, '
             f'which contains 90% of the points (points beyond it on the edge, triangles); grey line, x = y.\nRow 3: for '
             f'each cutoff, events with -log10 p.adjust above it in the within-species pair divided by those in '
             f'H2 vs M1 (same events).', ha='center', fontsize=6)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    for ext in ('pdf', 'png'):
        fig.savefig(OUT + f'Figure_padj_H2M1_vs_within.{ext}', dpi=300)
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(OUT + 'padj_pairs.tsv', sep='\t', index=False)
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()
