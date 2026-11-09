#!/usr/bin/env python3
"""
ISF figure: junction PSI of human spleen NK cells (GSE133383) against
  A  human blood NK, same dataset (GSE133383)
  B  human blood NK, another dataset (GSE115736, H2)
  C  mouse spleen NK (GSE233347; junctions mapped to the orthologous human coordinates)
for the junctions of the clusters with LeafCutter p.adjust <= 0.01 in GSE133383 spleen vs blood NK.
PSI per sample = junction reads / cluster reads (samples with >= 20 cluster reads); a point is the mean over the
samples of a group. Output: ensembl116/ISF/NK_spleen_PSI.pdf/.png/.tsv
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from effect_size_analysis import MIN_COVERAGE

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50/'
OUT = B + 'ISF/'
PADJ = 0.01
GROUPS = {
    'spleen': ('Human spleen NK (GSE133383)', 'NK_spleen_.*_GSE133383'),
    'blood': ('Human blood NK (GSE133383)', 'NK_blood_.*_GSE133383'),
    'blood_H2': ('Human blood NK (GSE115736)', r'NK_\d+_GSE115736'),
    'mouse_spleen': ('Mouse spleen NK (GSE233347)', 'NK_spleen_.*_GSE233347'),
}
matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 8, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})


def main():
    sig = pd.read_csv(ED + 'leafcutter/confounders/human_tissue/NK/leafcutter_ds_cluster_significance.txt', sep='\t')
    sig = set(sig.loc[(sig['status'] == 'Success') & (sig['p.adjust'] <= PADJ), 'cluster'])
    counts = pd.read_csv(ED + 'AS_clusters_value_confounders.txt', sep='\t', index_col=0)
    parts = counts.index.str.split(':')
    counts['cluster'] = [f'{p[0]}:{p[3]}' for p in parts]
    counts = counts[counts['cluster'].isin(sig)]
    psi = pd.DataFrame(index=counts.index)
    n = {}
    for g, (_, pattern) in GROUPS.items():
        samples = [c for c in counts.columns if pd.Series(c).str.fullmatch(pattern + '.*').iloc[0]]
        totals = counts.groupby('cluster')[samples].transform('sum')
        psi[g] = (counts[samples] / totals.where(totals >= MIN_COVERAGE)).mean(axis=1)
        n[g] = len(samples)
    psi.insert(0, 'cluster', counts['cluster'])
    psi.index.name = 'junction'
    os.makedirs(OUT, exist_ok=True)
    psi.round(4).to_csv(OUT + 'NK_spleen_PSI.tsv', sep='\t')

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw={'wspace': 0.4})
    for ax, g, letter in zip(axes, ('blood', 'blood_H2', 'mouse_spleen'), 'ABC'):
        d = psi[['spleen', g]].dropna()
        r = stats.pearsonr(d['spleen'], d[g])[0]
        ax.plot([0, 1], [0, 1], color='#9a9a9a', lw=0.8, ls='--', zorder=0)
        ax.scatter(d['spleen'], d[g], s=9, color='#2c5f8e', alpha=0.6, linewidths=0)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        ax.set_xticks([0, 0.5, 1])
        ax.set_yticks([0, 0.5, 1])
        ax.set_xlabel(f'PSI, {GROUPS["spleen"][0]}\n(n = {n["spleen"]})', fontsize=7)
        ax.set_ylabel(f'PSI, {GROUPS[g][0]}\n(n = {n[g]})', fontsize=7)
        ax.text(0.03, 0.97, f'r = {r:.2f}\n{len(d)} junctions\n{d.join(psi["cluster"])["cluster"].nunique()} clusters',
                transform=ax.transAxes, ha='left', va='top', fontsize=6.5,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.5), zorder=3)
        ax.text(-0.32, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight='bold')
        print(g, 'junctions', len(d), 'r', round(r, 3))
    fig.savefig(OUT + 'NK_spleen_PSI.pdf', bbox_inches='tight')
    fig.savefig(OUT + 'NK_spleen_PSI.png', dpi=300, bbox_inches='tight')
    print('clusters', len(sig), 'junctions', len(psi), 'samples', n)


if __name__ == '__main__':
    main()
