#!/usr/bin/env python3
"""
Candidate figure: sample-by-sample correlation matrix of junction PSI, all immune samples of the five datasets
(H1, H2, H3, M1, M2).

PSI of a junction in a sample = junction reads / reads of its cluster (the AS event), when the cluster has at
least 20 reads in that sample (as in the classification). Each cell = Pearson r over the orthologous junctions
with a PSI in both samples (pairwise-complete; the number of junctions is written to the TSV).
Ordering: the replicates of a dataset x cell type stay together. The groups (dataset x cell type) are ordered by
average-linkage hierarchical clustering of 1 - mean r between groups (dendrogram on top); the samples inside a
group are ordered by clustering their own sub-matrix.

Junction sets (argument; default: all three):
  all     every junction of every event, pairwise complete
  one     one junction per event (the first junction of the event in genomic order), pairwise complete
  common  every junction with a PSI in all 104 samples (the same junctions for every pair)

Output: ensembl116/revision_separation/candidate_figures/Figure_sample_correlation_<set>.(pdf|png),
        sample_psi_correlation_<set>.tsv (r), sample_psi_correlation_<set>_n.tsv (junctions per pair),
        sample_psi_correlation_summary.tsv (median r per pair of datasets and per pair type)
Usage: figure_sample_correlation_116.py [all|one|common ...]
"""
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import squareform

matplotlib.rcParams['font.family'] = 'Liberation Sans'
matplotlib.rcParams['pdf.fonttype'] = 42

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50_sep/'
OUT = B + 'revision_separation/candidate_figures/'
DATASETS = ['H1', 'H2', 'H3', 'M1', 'M2']
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_LABEL = {'CD4T': 'T4', 'CD8T': 'T8', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo'}
MIN_READS = 20
SPECIES_COLOR = {'Human': '#2c5f8e', 'Mouse': '#b4541a'}
DATASET_COLOR = {'H1': '#9ecae1', 'H2': '#4292c6', 'H3': '#08519c', 'M1': '#fdae6b', 'M2': '#d94801'}
TITLE = {'all': 'all junctions', 'one': 'one junction per event', 'common': 'junctions in all samples'}
CELL_COLOR = {'CD4T': '#1b9e77', 'CD8T': '#7570b3', 'NveB': '#e7298a', 'NK': '#66a61e', 'Mono': '#e6ab02'}


def sample_psi():
    design = pd.read_csv(ED + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(DATASETS)]
    counts = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', index_col=0)[design.index]
    cluster = counts.index.str.split(':').str[3]
    tot = counts.groupby(cluster).transform('sum')
    psi = counts / tot.where(tot >= MIN_READS)
    return psi.dropna(how='all'), design


def order_samples(r, design):
    design = design.assign(group=design['dataset'] + '_' + design['cell'])
    groups = sorted(design['group'].unique())
    members = {g: design.index[design['group'] == g].tolist() for g in groups}
    gr = pd.DataFrame([[r.loc[members[a], members[b]].to_numpy().mean() for b in groups] for a in groups],
                      index=groups, columns=groups)
    d = 1 - gr.to_numpy()
    np.fill_diagonal(d, 0)
    link = linkage(squareform((d + d.T) / 2, checks=False), 'average')
    order = []
    for g in [groups[i] for i in leaves_list(link)]:
        s = members[g]
        if len(s) > 2:
            sub = 1 - r.loc[s, s].to_numpy()
            np.fill_diagonal(sub, 0)
            s = [s[i] for i in leaves_list(linkage(squareform((sub + sub.T) / 2, checks=False), 'average'))]
        order += s
    return order, link, groups, members


def plot(r, design, order, link, groups, members, mode, n_junctions):
    n = len(order)
    fig = plt.figure(figsize=(9.2, 10.4))
    gs = fig.add_gridspec(5, 3, height_ratios=[1.3, 0.16, 0.16, 0.16, 8], width_ratios=[8, 0.25, 0.25],
                          hspace=0.03, wspace=0.05, left=0.1, right=0.9, top=0.97, bottom=0.085)
    # dendrogram of the groups, leaves placed at the centre of each group block
    axd = fig.add_subplot(gs[0, 0])
    dn = dendrogram(link, no_plot=True, labels=groups)
    start, centre = 0, {}
    for g in dn['ivl']:
        centre[g] = start + len(members[g]) / 2
        start += len(members[g])
    leaf_x = [centre[g] for g in dn['ivl']]
    xmap = lambda x: np.interp(x, 5 + 10 * np.arange(len(groups)), leaf_x)
    for xs, ys in zip(dn['icoord'], dn['dcoord']):
        axd.plot([xmap(x) for x in xs], ys, color='#555555', lw=0.8)
    axd.set_xlim(0, n)
    axd.axis('off')
    axd.set_title(f'Correlation of junction PSI between samples (Pearson r) - {TITLE[mode]} ({n_junctions:,} junctions)',
                  fontsize=10, loc='left')
    # annotation bars
    ann = [('Species', design.loc[order, 'dataset'].str[0].map({'H': 'Human', 'M': 'Mouse'}).map(SPECIES_COLOR)),
           ('Dataset', design.loc[order, 'dataset'].map(DATASET_COLOR)),
           ('Cell type', design.loc[order, 'cell'].map(CELL_COLOR))]
    for k, (name, cols) in enumerate(ann):
        ax = fig.add_subplot(gs[1 + k, 0])
        ax.imshow([matplotlib.colors.to_rgba_array(cols.tolist())], aspect='auto', extent=(0, n, 0, 1))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(name, rotation=0, ha='right', va='center', fontsize=7)
        for s in ax.spines.values():
            s.set_visible(False)
    # matrix
    ax = fig.add_subplot(gs[4, 0])
    m = r.loc[order, order].to_numpy()
    vmin = np.floor(np.nanmin(m) * 20) / 20
    im = ax.imshow(m, cmap='viridis', vmin=vmin, vmax=1, extent=(0, n, n, 0), interpolation='nearest')
    edges = np.cumsum([len(members[g]) for g in dn['ivl']])
    for e in edges[:-1]:
        ax.axhline(e, color='white', lw=0.6)
        ax.axvline(e, color='white', lw=0.6)
    ticks = [centre[g] for g in dn['ivl']]
    labels = [f"{g.split('_')[0]} {CELL_LABEL[g.split('_')[1]]}" for g in dn['ivl']]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.tick_params(length=0)
    # colour bar
    cax = fig.add_subplot(gs[4, 2])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label('Pearson r', fontsize=8)
    cb.ax.tick_params(labelsize=7)
    handles = ([Patch(color=c, label=k) for k, c in SPECIES_COLOR.items()] +
               [Patch(color=c, label=k) for k, c in DATASET_COLOR.items()] +
               [Patch(color=CELL_COLOR[c], label=CELL_LABEL[c]) for c in CELLS])
    fig.legend(handles=handles, loc='lower center', ncol=12, fontsize=7, frameon=False,
               bbox_to_anchor=(0.5, 0.005), handlelength=1, columnspacing=0.9)
    for ext in ('pdf', 'png'):
        fig.savefig(OUT + f'Figure_sample_correlation_{mode}.{ext}', dpi=300)
    plt.close(fig)


def select(psi, mode):
    if mode == 'one':
        clu = psi.index.str.split(':').str[3]
        start = psi.index.str.split(':').str[1].astype(int)
        first = pd.Series(start, index=psi.index).groupby(clu).idxmin()
        return psi.loc[first.to_numpy()]
    if mode == 'common':
        return psi.dropna()
    return psi


def pair_summary(r, design, mode):
    ds = design.loc[r.index, 'dataset'].to_numpy()
    cell = design.loc[r.index, 'cell'].to_numpy()
    i, j = np.triu_indices(len(r), 1)
    v = r.to_numpy()[i, j]
    sp_i, sp_j = np.array([d[0] for d in ds[i]]), np.array([d[0] for d in ds[j]])
    kind = np.where(sp_i != sp_j, 'human vs mouse',
                    np.where(ds[i] != ds[j], 'same species, other dataset',
                             np.where(cell[i] == cell[j], 'same dataset, same cell type', 'same dataset, other cell type')))
    pair = np.array(['-'.join(sorted(x)) for x in zip(ds[i], ds[j])])
    rows = [pd.Series(v).groupby(kind).median(), pd.Series(v).groupby(pair).median()]
    return pd.concat(rows).rename(mode)


def run(psi_all, design, mode):
    psi = select(psi_all, mode)
    r = psi.corr(method='pearson', min_periods=100)
    present = psi.notna().astype(float)
    nj = present.T @ present
    order, link, groups, members = order_samples(r, design)
    r.loc[order, order].to_csv(OUT + f'sample_psi_correlation_{mode}.tsv', sep='\t', float_format='%.4f')
    nj.loc[order, order].astype(int).to_csv(OUT + f'sample_psi_correlation_{mode}_n.tsv', sep='\t')
    plot(r, design, order, link, groups, members, mode, len(psi))
    off = nj.to_numpy()[~np.eye(len(nj), dtype=bool)]
    print(f'{mode}: {len(psi)} junctions; per pair median {np.median(off):.0f}, min {off.min():.0f}; '
          f'r {np.nanmin(r.to_numpy()):.3f}-1')
    print('  group order:', dendrogram(link, no_plot=True, labels=groups)['ivl'])
    return pair_summary(r, design, mode)


def main():
    os.makedirs(OUT, exist_ok=True)
    modes = sys.argv[1:] or ['all', 'one', 'common']
    psi, design = sample_psi()
    summary = pd.concat([run(psi, design, m) for m in modes], axis=1)
    summary.index.name = 'pairs'
    summary.round(3).to_csv(OUT + 'sample_psi_correlation_summary.tsv', sep='\t')
    print(summary.round(3))


if __name__ == '__main__':
    main()
