#!/usr/bin/env python3
"""
ISF figure: junction dPSI over all clusters (orthologous junctions, human coordinates).
  x (both panels)  human spleen NK GSE133383 vs mouse spleen NK GSE233347
  y, panel A       human spleen NK GSE133383 vs human blood NK GSE133383
  y, panel B       human blood NK GSE133383 vs human blood NK GSE115736
dPSI of "P vs Q" = PSI(Q) - PSI(P); PSI as in isf_nk_spleen_psi_116.py (mean of per-sample PSI, samples with
>= 20 cluster reads, or --min-reads). A junction is shown when all groups of its panel have PSI.
Output: ensembl116/ISF/NK_dPSI[_min<N>reads].pdf/.png/.tsv
"""
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

from effect_size_analysis import MIN_COVERAGE
from isf_nk_spleen_psi_116 import ED, GROUPS, OUT

COMPARISONS = {
    'spleen_vs_mouse_spleen': ('spleen', 'mouse_spleen'),
    'spleen_vs_blood': ('spleen', 'blood'),
    'blood_vs_blood_H2': ('blood', 'blood_H2'),
}
SHORT = {'spleen': 'human spleen (GSE133383)', 'blood': 'human blood (GSE133383)',
         'blood_H2': 'human blood (GSE115736)', 'mouse_spleen': 'mouse spleen (GSE233347)'}


def group_psi(min_reads=MIN_COVERAGE):
    counts = pd.read_csv(ED + 'AS_clusters_value_confounders.txt', sep='\t', index_col=0)
    cluster = [f'{p[0]}:{p[3]}' for p in counts.index.str.split(':')]
    psi = pd.DataFrame(index=counts.index)
    for g, (_, pattern) in GROUPS.items():
        samples = counts.columns[counts.columns.str.fullmatch(pattern + '.*')].tolist()
        totals = counts[samples].groupby(cluster).transform('sum')
        psi[g] = (counts[samples] / totals.where(totals >= min_reads)).mean(axis=1)
    psi.insert(0, 'cluster', cluster)
    psi.index.name = 'junction'
    return psi


def main(min_reads):
    psi = group_psi(min_reads)
    tag = '' if min_reads == MIN_COVERAGE else f'_min{min_reads}reads'
    d = psi[['cluster']].copy()
    for name, (p, q) in COMPARISONS.items():
        d[f'dPSI_{name}'] = psi[q] - psi[p]
    d.round(4).to_csv(OUT + f'NK_dPSI{tag}.tsv', sep='\t')

    x = 'dPSI_spleen_vs_mouse_spleen'
    fig, axes = plt.subplots(1, 2, figsize=(5.4, 2.8), gridspec_kw={'wspace': 0.45})
    for ax, y, letter in zip(axes, ('dPSI_spleen_vs_blood', 'dPSI_blood_vs_blood_H2'), 'AB'):
        sub = d[[x, y, 'cluster']].dropna()
        r = stats.pearsonr(sub[x], sub[y])[0]
        for v in (ax.axhline, ax.axvline):
            v(0, color='#9a9a9a', lw=0.6, zorder=0)
        ax.scatter(sub[x], sub[y], s=3, color='#2c5f8e', alpha=0.2, linewidths=0, rasterized=True)
        ax.set_xlim(-1.02, 1.02)
        ax.set_ylim(-1.02, 1.02)
        ax.set_aspect('equal')
        ticks = [-1, -0.5, 0, 0.5, 1]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        p, q = COMPARISONS[x[5:]]
        ax.set_xlabel(f'ΔPSI, {SHORT[q]}\n− {SHORT[p]}', fontsize=7)
        p, q = COMPARISONS[y[5:]]
        ax.set_ylabel(f'ΔPSI, {SHORT[q]}\n− {SHORT[p]}', fontsize=7)
        ax.text(0.03, 0.97, f'r = {r:.2f}\n{len(sub):,} junctions\n{sub["cluster"].nunique():,} clusters',
                transform=ax.transAxes, ha='left', va='top', fontsize=6.5,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.5), zorder=3)
        ax.text(-0.3, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight='bold')
        print(y, 'junctions', len(sub), 'clusters', sub['cluster'].nunique(), 'r', round(r, 3),
              'sd x', round(sub[x].std(), 3), 'sd y', round(sub[y].std(), 3))
    fig.savefig(OUT + f'NK_dPSI{tag}.pdf', bbox_inches='tight', dpi=300)
    fig.savefig(OUT + f'NK_dPSI{tag}.png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 8, 'pdf.fonttype': 42,
                                'axes.spines.top': False, 'axes.spines.right': False})
    ap = argparse.ArgumentParser()
    ap.add_argument('--min-reads', type=int, default=MIN_COVERAGE,
                    help='minimum reads of the cluster in a sample for its PSI to be used (default 20)')
    main(ap.parse_args().min_reads)
