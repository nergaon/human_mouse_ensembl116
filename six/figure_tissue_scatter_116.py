#!/usr/bin/env python3
"""
Candidate figure: tissue controls. One panel per tissue comparison (same cell type, same dataset, two tissues);
each point = one event testable in the comparison, at the junction with the largest |dPSI| (the max|dPSI| junction);
x / y = mean PSI over the replicates of tissue 1 / tissue 2 (samples with at least 20 JSRs in the event, as for
max|dPSI|). Events differentially spliced in the control (LeafCutter p.adjust <= 0.05 and max|dPSI| >= 0.1) are
highlighted.
Usage: figure_tissue_scatter_116.py [events root]   (default ensembl116/events; the backup
       ensembl116/pre_liftover_20261118/events can be given)
Output: revision_6ds/figures/Figure_tissue_controls_scatter.pdf / .tif and .tsv (the plotted values)
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from figures_main_116 import CELL_LABEL, OUT, save

ROOT = (sys.argv[1].rstrip('/') + '/') if len(sys.argv) > 1 else \
    '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/'
TISSUE = [('human_tissue', 'Human NK cells', 'GSE133383'),
          ('mouse_tissue', 'Mouse NK cells', 'GSE233347'),
          ('mouse_tissue_m1', 'Mouse monocytes', 'GSE116177'),
          ('mouse_tissue_mono', 'Mouse monocytes', 'PRJNA797236')]
POINT, DS_COLOR, INK, MUTED = '#9aa5b1', '#c0392b', '#1f2328', '#6b7280'


def main():
    j = pd.read_csv(ROOT + 'id50_6ds/effects_confounders/junction_effects.tsv', sep='\t',
                    usecols=['comparison', 'cell', 'group1', 'group2', 'cluster', 'psi1', 'psi2', 'dpsi', 'testable'])
    ev = pd.read_csv(ROOT + 'id50_6ds_anova/effects_confounders/event_effects.tsv', sep='\t',
                     usecols=['comparison', 'cell', 'cluster', 'testable', 'status_m0.1'])
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 6.4), gridspec_kw={'hspace': 0.42, 'wspace': 0.32})
    rows = []
    for ax, (comp, label, ds) in zip(axes.flat, TISSUE):
        e = ev[(ev['comparison'] == comp) & ev['testable']]
        x = j[(j['comparison'] == comp) & j['cluster'].isin(e['cluster']) & j['dpsi'].notna()]
        x = x.loc[x['dpsi'].abs().groupby(x['cluster']).idxmax()]          # one junction per event
        x = x.merge(e[['cluster', 'status_m0.1']], on='cluster')
        ds_mask = x['status_m0.1'] == 'DS'
        g1, g2 = x['group1'].iloc[0].replace('_', ' '), x['group2'].iloc[0].replace('_', ' ')
        ax.plot([0, 1], [0, 1], color=MUTED, lw=0.6, ls='--', zorder=1)
        ax.scatter(x.loc[~ds_mask, 'psi1'], x.loc[~ds_mask, 'psi2'], s=4, color=POINT, alpha=0.5, lw=0, zorder=2,
                   rasterized=True)
        ax.scatter(x.loc[ds_mask, 'psi1'], x.loc[ds_mask, 'psi2'], s=10, color=DS_COLOR, lw=0.4, edgecolor='white',
                   zorder=3)
        r = np.corrcoef(x['psi1'], x['psi2'])[0, 1]
        ax.text(0.03, 0.97, f'n = {len(x):,} events\nr = {r:.3f}\nDS: {int(ds_mask.sum())}', transform=ax.transAxes,
                va='top', fontsize=6.5, color=INK)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        ax.set_xlabel(f'PSI, {g1} (mean of replicates)', fontsize=7)
        ax.set_ylabel(f'PSI, {g2} (mean of replicates)', fontsize=7)
        ax.set_title(f'{label}: {g1} vs {g2} ({ds})', fontsize=7)
        rows.append(x.assign(control=comp, dataset=ds))
    fig.legend([Line2D([], [], marker='o', ls='', color=POINT, ms=3),
                Line2D([], [], marker='o', ls='', color=DS_COLOR, ms=4)],
               ['Event (junction with the largest |ΔPSI|)', 'Differentially spliced in the control'],
               loc='upper center', ncol=2, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.975))
    save(fig, 'Figure_tissue_controls_scatter')
    pd.concat(rows)[['control', 'dataset', 'cell', 'cluster', 'group1', 'group2', 'psi1', 'psi2', 'dpsi',
                     'status_m0.1']].to_csv(OUT + 'Figure_tissue_controls_scatter.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
