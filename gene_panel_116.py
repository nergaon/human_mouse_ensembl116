#!/usr/bin/env python3
"""
Panel B of the gene figures: junction spanning read (JSR) counts (top) and PSI (bottom) per sample,
for one AS event, across the immune datasets and fibroblasts. Layout follows genes_figure.py
(ensembl115): cell-type blocks T4, T8, B, NK, Mo, then fibroblasts; datasets within a block in the
order M2, M1, H1, H2, H3 (fibroblasts: M, H).

Usage: gene_panel_116.py <cluster id, e.g. clu_52221> <output prefix> [--title USP16]
                         [--colors '#1f4e79,#9dc3e6']   (one per junction, in genomic order)
"""
import argparse
import re

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

matplotlib.rcParams['font.family'] = 'Liberation Sans'
matplotlib.rcParams['pdf.fonttype'] = 42

EVENTS = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50/'
CELL_ORDER = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Fibro']
CELL_DISPLAY = {'CD4T': 'T4', 'CD8T': 'T8', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo', 'Fibro': 'Fibroblasts'}
DATASET_ORDER = ['M2', 'M1', 'H1', 'H2', 'H3', 'FM', 'FH']
DATASET_DISPLAY = {'FM': 'M', 'FH': 'H'}


def load(cluster):
    counts = pd.read_csv(EVENTS + 'AS_clusters_value.txt', sep='\t', index_col=0)
    rows = counts[counts.index.str.endswith(':' + cluster)]
    design = pd.read_csv(EVENTS + 'leafcutter/sample_design.tsv', sep='\t')
    design = design[design['sample'].isin(rows.columns)]
    return rows, design


def natural_key(text):
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r'(\d+)', text)]


def positions_for(samples, design, bar_w, dataset_gap=0.5, cell_gap=0.65, fibro_gap=2.5):
    cell = design.set_index('sample')['cell']
    ds = design.set_index('sample')['dataset']
    pos = np.zeros(len(samples))
    for i in range(1, len(samples)):
        pos[i] = pos[i - 1] + bar_w
        if ds[samples[i]] != ds[samples[i - 1]]:
            pos[i] += dataset_gap
        if cell[samples[i]] != cell[samples[i - 1]]:
            pos[i] += fibro_gap if 'Fibro' in (cell[samples[i]], cell[samples[i - 1]]) else cell_gap
    return pos


def blocks(samples, key):
    out, start = [], 0
    for i in range(1, len(samples) + 1):
        if i == len(samples) or key[samples[i]] != key[samples[start]]:
            out.append((key[samples[start]], start, i - 1))
            start = i
    return out


def plot(cluster, out_prefix, title=None, colors=None):
    rows, design = load(cluster)
    design = design.assign(c=design['cell'].map(CELL_ORDER.index), d=design['dataset'].map(DATASET_ORDER.index),
                           n=design['sample'].map(lambda s: tuple(map(str, natural_key(s)))))
    design = design.sort_values(['c', 'd', 'n'])
    samples = design['sample'].tolist()
    cell = design.set_index('sample')['cell']
    ds = design.set_index('sample')['dataset']
    cnt = rows[samples].to_numpy(float)
    psi = (rows[samples] / rows[samples].sum(axis=0).replace(0, np.nan)).fillna(0).to_numpy(float)
    labels = [f"{j.split(':')[0]}:{j.split(':')[1]}-{j.split(':')[2]}" for j in rows.index]
    if colors is None:
        cmap = plt.get_cmap('tab20')
        colors = [cmap(i % cmap.N) for i in range(len(rows))]

    bar_w = 0.19
    pos = positions_for(samples, design, bar_w)
    cell_blocks = blocks(samples, cell)
    plt.close('all')
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(11.69, 4.13), gridspec_kw={'height_ratios': [1, 1]})
    b_cnt, b_psi = np.zeros(len(samples)), np.zeros(len(samples))
    for i in range(len(rows)):
        top.bar(pos, cnt[i], width=bar_w, bottom=b_cnt, color=colors[i], edgecolor='none')
        bottom.bar(pos, psi[i], width=bar_w, bottom=b_psi, color=colors[i], edgecolor='none')
        b_cnt += cnt[i]
        b_psi += psi[i]

    left = right = bar_w * 0.55
    if len(cell_blocks) > 1:
        left = max(left, (pos[cell_blocks[1][1]] - pos[cell_blocks[0][2]]) / 2)
        right = max(right, (pos[cell_blocks[-1][1]] - pos[cell_blocks[-2][2]]) / 2)
    for ax in (top, bottom):
        ax.set_xlim(pos.min() - left, pos.max() + right)

    top.set_xticks([pos[s:e + 1].mean() for _, s, e in cell_blocks])
    top.set_xticklabels([CELL_DISPLAY[c] for c, _, _ in cell_blocks], fontsize=10)
    top.xaxis.tick_top()
    top.tick_params(axis='x', length=0, labeltop=True, labelbottom=False)
    top.set_ylabel('JSR counts')

    # dataset labels: contiguous dataset runs within each cell block
    ticks, tick_labels = [], []
    for c, s, e in cell_blocks:
        sub = samples[s:e + 1]
        start = 0
        for i in range(1, len(sub) + 1):
            if i == len(sub) or ds[sub[i]] != ds[sub[start]]:
                ticks.append(pos[s + start:s + i].mean())
                tick_labels.append(DATASET_DISPLAY.get(ds[sub[start]], ds[sub[start]]))
                start = i
    bottom.set_xticks(ticks)
    bottom.set_xticklabels(tick_labels, fontsize=8)
    bottom.tick_params(axis='x', length=0)
    bottom.set_ylabel('PSI')
    bottom.set_ylim(0, 1)

    for i in range(1, len(samples)):
        if ds[samples[i]] != ds[samples[i - 1]] and cell[samples[i]] == cell[samples[i - 1]]:
            x = (pos[i - 1] + pos[i]) / 2
            for ax in (top, bottom):
                ax.axvline(x, color='#aaaaaa', linewidth=0.7, linestyle='--', zorder=0)
    for (c, s, e), (c2, _, _) in zip(cell_blocks[:-1], cell_blocks[1:]):
        x = (pos[e] + pos[e + 1]) / 2
        fibro = 'Fibro' in (c, c2)
        for ax in (top, bottom):
            if fibro:
                ax.axvspan(x - bar_w * 1.2, x + bar_w * 1.2, color='#f0f0f0', alpha=0.9, zorder=0)
            ax.axvline(x, color='#4d4d4d' if fibro else '#999999', linewidth=1.8 if fibro else 0.8)

    fig.suptitle(title or cluster, fontsize=10)
    fig.legend([Patch(facecolor=colors[i], edgecolor='none') for i in range(len(rows))], labels, fontsize=7,
               loc='upper center', bbox_to_anchor=(0.5, -0.03), ncol=4, frameon=False)
    plt.subplots_adjust(hspace=0.30)
    fig.savefig(out_prefix + '.pdf', bbox_inches='tight', dpi=300)
    fig.savefig(out_prefix + '.png', bbox_inches='tight', dpi=300)
    plt.close('all')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cluster')
    ap.add_argument('out_prefix')
    ap.add_argument('--title')
    ap.add_argument('--colors', help='comma separated, one per junction in genomic order')
    a = ap.parse_args()
    plot(a.cluster, a.out_prefix, a.title, a.colors.split(',') if a.colors else None)


if __name__ == '__main__':
    main()
