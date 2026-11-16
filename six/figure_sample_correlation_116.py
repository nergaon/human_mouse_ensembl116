#!/usr/bin/env python3
"""
Candidate figure: sample-by-sample correlation matrix of junction PSI, all immune samples of the five datasets
(H1, H2, H3, M1, M2).

PSI of a junction in a sample = junction reads / reads of its cluster (the AS event), when the cluster has at
least 20 reads in that sample (as in the classification). Each cell = Pearson r over the orthologous junctions
with a PSI in both samples (pairwise-complete; the number of junctions is written to the TSV).
Ordering (sorted, no clustering): dataset (H1, H2, H3, M1, M2), then cell type in a fixed order (T4, T8, B, NK, Mo),
then sample name, so that a cell-type effect shows as the same pattern in every dataset-by-dataset block.
Cell-type effect: for pairs of samples from different datasets of the same species, and from human vs mouse, the
mean r of same-cell-type pairs minus the mean r of other-cell-type pairs, averaged over the dataset pairs; P from
1,000 permutations of the cell-type labels between the dataset x cell type groups of each dataset (replicates move
together). Within a dataset, same cell type = replicates, so no permutation P is given there.

Junction sets (argument; default: all four):
  all     every junction of every event, pairwise complete
  one     one junction per event (the first junction of the event in genomic order), pairwise complete
  common  every junction with a PSI in all 104 samples (the same junctions for every pair)
  ds      events differentially spliced (DS) in at least one cell type, one junction per event: the junction with
          the highest |dPSI| (mean PSI of H1-H3 minus mean PSI of M1-M2) in a cell type where the event is DS with
          high confidence; if none, in a cell type where it is DS (any confidence). The same junction is used in all
          cell types, pairwise complete

Output: ensembl116/revision_6ds/candidate_figures/Figure_sample_correlation_<set>.(pdf|png),
        sample_psi_correlation_<set>.tsv (r), sample_psi_correlation_<set>_n.tsv (junctions per pair),
        sample_psi_correlation_summary.tsv (median r per pair of datasets and per pair type)
Usage: figure_sample_correlation_116.py [all|one|common|ds ...]
"""
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

matplotlib.rcParams['font.family'] = 'Liberation Sans'
matplotlib.rcParams['pdf.fonttype'] = 42

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50_6ds_anova/'
OUT = B + 'revision_6ds/candidate_figures/'
DATASETS = ['H1', 'H2', 'H3', 'M1', 'M2']
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_LABEL = {'CD4T': 'T4', 'CD8T': 'T8', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo'}
MIN_READS = 20
SPECIES_COLOR = {'Human': '#2c5f8e', 'Mouse': '#b4541a'}
DATASET_COLOR = {'H1': '#9ecae1', 'H2': '#4292c6', 'H3': '#08519c', 'M1': '#fdae6b', 'M2': '#d94801'}
TITLE = {'all': 'all junctions', 'one': 'one junction per event', 'common': 'junctions in all samples',
         'ds': 'one junction per DS event'}
CELL_COLOR = {'CD4T': '#1b9e77', 'CD8T': '#7570b3', 'NveB': '#e7298a', 'NK': '#66a61e', 'Mono': '#e6ab02'}


def sample_psi():
    design = pd.read_csv(ED + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(DATASETS)]
    counts = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', index_col=0)[design.index]
    cluster = counts.index.str.split(':').str[3]
    tot = counts.groupby(cluster).transform('sum')
    psi = counts / tot.where(tot >= MIN_READS)
    return psi.dropna(how='all'), design


def order_samples(design):
    """Fixed order (no clustering): dataset (H1, H2, H3, M1, M2), then cell type (CELLS), then sample name."""
    d = design.assign(ds_i=design['dataset'].map(DATASETS.index), cell_i=design['cell'].map(CELLS.index))
    d = d.sort_values(['ds_i', 'cell_i']).rename_axis('sample').reset_index()
    d = d.sort_values(['ds_i', 'cell_i', 'sample'])
    groups = list(dict.fromkeys(d['dataset'] + '_' + d['cell']))
    members = {g: d.loc[d['dataset'] + '_' + d['cell'] == g, 'sample'].tolist() for g in groups}
    return d['sample'].tolist(), groups, members


def plot(r, design, order, groups, members, mode, n_junctions):
    n = len(order)
    fig = plt.figure(figsize=(9.2, 9.4))
    gs = fig.add_gridspec(4, 3, height_ratios=[0.16, 0.16, 0.16, 8], width_ratios=[8, 0.25, 0.25],
                          hspace=0.03, wspace=0.05, left=0.1, right=0.9, top=0.95, bottom=0.095)
    # annotation bars
    ann = [('Species', design.loc[order, 'dataset'].str[0].map({'H': 'Human', 'M': 'Mouse'}).map(SPECIES_COLOR)),
           ('Dataset', design.loc[order, 'dataset'].map(DATASET_COLOR)),
           ('Cell type', design.loc[order, 'cell'].map(CELL_COLOR))]
    for k, (name, cols) in enumerate(ann):
        ax = fig.add_subplot(gs[k, 0])
        ax.imshow([matplotlib.colors.to_rgba_array(cols.tolist())], aspect='auto', extent=(0, n, 0, 1))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(name, rotation=0, ha='right', va='center', fontsize=7)
        for sp in ax.spines.values():
            sp.set_visible(False)
        if k == 0:
            ax.set_title(f'Correlation of junction PSI between samples (Pearson r) - {TITLE[mode]} '
                         f'({n_junctions:,} junctions)', fontsize=10, loc='left', pad=8)
    # matrix
    ax = fig.add_subplot(gs[3, 0])
    m = r.loc[order, order].to_numpy()
    vmin = np.floor(np.nanmin(m) * 20) / 20
    im = ax.imshow(m, cmap='viridis', vmin=vmin, vmax=1, extent=(0, n, n, 0), interpolation='nearest')
    sizes = [len(members[g]) for g in groups]
    edges = np.cumsum(sizes)
    for e, g in zip(edges[:-1], groups):
        new_ds = g.split('_')[0] != groups[groups.index(g) + 1].split('_')[0]
        lw, col = (1.6, 'white') if new_ds else (0.4, '#dddddd')
        ax.axhline(e, color=col, lw=lw)
        ax.axvline(e, color=col, lw=lw)
    ticks = edges - np.array(sizes) / 2
    labels = [f"{g.split('_')[0]} {CELL_LABEL[g.split('_')[1]]}" for g in groups]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.tick_params(length=0)
    # colour bar
    cax = fig.add_subplot(gs[3, 2])
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


def ds_junctions():
    """One junction per event DS in >= 1 cell type (see docstring)."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import classify_separation_116 as cs
    st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t',
                     usecols=['cell', 'cluster', 'status_high', 'status_high_low'])
    st = st[st['status_high_low'] == 'divergent'].assign(high=lambda d: d['status_high'] == 'divergent')
    counts = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(ED + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    rows = []
    for cell, d in st.groupby('cell'):
        psi, clu = cs.dataset_psi(counts, design, cell)
        dpsi = (psi[list(cs.HUMAN)].mean(axis=1) - psi[['M1', 'M2']].mean(axis=1)).abs()
        dpsi = pd.DataFrame({'junction': psi.index, 'cluster': clu, 'abs_dpsi': dpsi.to_numpy()})
        rows.append(dpsi.merge(d, on='cluster'))
    j = pd.concat(rows, ignore_index=True)
    j = j.sort_values(['cluster', 'high', 'abs_dpsi'], ascending=[True, False, False]).drop_duplicates('cluster')
    j[['cluster', 'junction', 'cell', 'high', 'abs_dpsi']].to_csv(OUT + 'sample_psi_correlation_ds_junctions.tsv',
                                                                sep='\t', index=False, float_format='%.4f')
    print(f"ds: {len(j)} DS events ({j['high'].sum()} chosen from a high-confidence cell type)")
    return j['junction'].to_numpy()


def select(psi, mode):
    if mode == 'ds':
        return psi.loc[psi.index.intersection(ds_junctions())]
    if mode == 'one':
        clu = psi.index.str.split(':').str[3]
        start = psi.index.str.split(':').str[1].astype(int)
        first = pd.Series(start, index=psi.index).groupby(clu).idxmin()
        return psi.loc[first.to_numpy()]
    if mode == 'common':
        return psi.dropna()
    return psi


def pair_kinds(r, design):
    ds = design.loc[r.index, 'dataset'].to_numpy()
    cell = design.loc[r.index, 'cell'].to_numpy()
    i, j = np.triu_indices(len(r), 1)
    sp_i, sp_j = np.array([d[0] for d in ds[i]]), np.array([d[0] for d in ds[j]])
    kind = np.where(sp_i != sp_j, 'human vs mouse',
                    np.where(ds[i] != ds[j], 'same species, other dataset', 'same dataset'))
    pair = np.array(['-'.join(sorted(x)) for x in zip(ds[i], ds[j])])
    return i, j, kind, pair, ds, cell


def cell_effect(v, i, j, kind, pair, ds, cell, k):
    """Mean over the dataset pairs of kind k of (mean r same cell type - mean r other cell type)."""
    sel = kind == k
    same = cell[i] == cell[j]
    df = pd.DataFrame({'pair': pair[sel], 'same': same[sel], 'r': v[sel]}).dropna()
    m = df.groupby(['pair', 'same'])['r'].mean().unstack()
    return (m[True] - m[False]).mean()


def pair_summary(r, design, mode, n_perm=1000, seed=1):
    i, j, kind, pair, ds, cell = pair_kinds(r, design)
    v = r.to_numpy()[i, j]
    same = np.where(cell[i] == cell[j], 'same cell type', 'other cell type')
    rows = {}
    for k in ['same dataset', 'same species, other dataset', 'human vs mouse']:
        for c in ['same cell type', 'other cell type']:
            rows[f'{k}, {c}'] = np.nanmedian(v[(kind == k) & (same == c)])
    rows['human vs mouse'] = np.nanmedian(v[kind == 'human vs mouse'])
    rows['same species, other dataset'] = np.nanmedian(v[kind == 'same species, other dataset'])
    for p in sorted(set(pair)):
        rows[p] = np.nanmedian(v[pair == p])
    # cell-type effect and permutation P: cell-type labels permuted between the dataset x cell type groups of
    # each dataset (replicates move together)
    rng = np.random.default_rng(seed)
    grp = pd.Series(ds + '_' + cell)
    rows['same dataset: cell-type effect (r same - r other)'] = cell_effect(v, i, j, kind, pair, ds, cell,
                                                                            'same dataset')
    for k in ['same species, other dataset', 'human vs mouse']:
        obs = cell_effect(v, i, j, kind, pair, ds, cell, k)
        null = []
        for _ in range(n_perm):
            relabel = {}
            for d in DATASETS:
                cs_ = [c for c in CELLS if f'{d}_{c}' in set(grp)]
                for a_, b_ in zip(cs_, rng.permutation(cs_)):
                    relabel[f'{d}_{a_}'] = b_
            pc = grp.map(relabel).to_numpy()
            null.append(cell_effect(v, i, j, kind, pair, ds, pc, k))
        rows[f'{k}: cell-type effect (r same - r other)'] = obs
        rows[f'{k}: cell-type effect P'] = (1 + np.sum(np.array(null) >= obs)) / (1 + n_perm)
    return pd.Series(rows).rename(mode)


def run(psi_all, design, mode):
    psi = select(psi_all, mode)
    r = psi.corr(method='pearson', min_periods=100)
    present = psi.notna().astype(float)
    nj = present.T @ present
    order, groups, members = order_samples(design.loc[r.index])
    r.loc[order, order].to_csv(OUT + f'sample_psi_correlation_{mode}.tsv', sep='\t', float_format='%.4f')
    nj.loc[order, order].astype(int).to_csv(OUT + f'sample_psi_correlation_{mode}_n.tsv', sep='\t')
    plot(r, design, order, groups, members, mode, len(psi))
    off = nj.to_numpy()[~np.eye(len(nj), dtype=bool)]
    print(f'{mode}: {len(psi)} junctions; per pair median {np.median(off):.0f}, min {off.min():.0f}; '
          f'r {np.nanmin(r.to_numpy()):.3f}-1')
    return pair_summary(r, design, mode)


def main():
    os.makedirs(OUT, exist_ok=True)
    modes = sys.argv[1:] or ['all', 'one', 'common', 'ds']
    psi, design = sample_psi()
    summary = pd.concat([run(psi, design, m) for m in modes], axis=1)
    summary.index.name = 'pairs'
    summary.round(3).to_csv(OUT + 'sample_psi_correlation_summary.tsv', sep='\t')
    print(summary.round(3))


if __name__ == '__main__':
    main()
