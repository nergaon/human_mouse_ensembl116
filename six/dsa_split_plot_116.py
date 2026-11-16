#!/usr/bin/env python3
"""Plots of DSA options 1 and 2 (option 1: Wilcoxon significant at FDR 5% and 10%, with and without the
requirement, for every split; option 2: number of stage-1 events separated by each split of the datasets, per cell type, for gap > 0.1
and gap > 0.2, with and without the requirement; the human-mouse split is highlighted. The plot is saved as PNG/PDF
and inserted into the option2_splits sheet of DSA_tables_<name>.xlsx.
Usage: dsa_split_plot_116.py [name]   (default 5datasets; run after dsa_analysis_116.py)"""
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from openpyxl.drawing.image import Image

D = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_6ds/dsa/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_NAME = {'CD4T': 'CD4 T', 'CD8T': 'CD8 T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Monocytes'}
COLS = ['gap > 0.1', 'gap > 0.1 + requirement', 'gap > 0.2', 'gap > 0.2 + requirement']
matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})


def plot(o, cols, title, png, rowlab=None):
    fig, axes = plt.subplots(len(cols), len(CELLS), figsize=(11, 2.1 * len(cols) + 0.5), sharey='row')
    for r, (col, sel) in enumerate(cols):
        for c, cell in enumerate(CELLS):
            ax = axes[r, c]
            base = o[o['cell'] == cell]
            order = base[sel(base)].sort_values(['human vs mouse', cols[0][0]], ascending=[False, False])['split']
            d = base[sel(base)].set_index('split').loc[order].reset_index()
            colors = ['#e8762c' if hm else '#9e9e9e' for hm in d['human vs mouse']]
            ax.bar(range(len(d)), d[col], color=colors, width=0.8)
            mixed = d.loc[~d['human vs mouse'], col]
            hm = d.loc[d['human vs mouse'], col].iloc[0]
            ax.axhline(mixed.mean(), color='black', ls='--', lw=0.7)
            ax.text(0, hm, f'{hm}', ha='center', va='bottom', fontsize=6)
            ax.set_xticks(range(len(d)))
            if r == len(cols) - 1:
                ax.set_xticklabels(d['split'].str.replace(' vs ', '\n'), rotation=90, fontsize=4.5)
            else:
                ax.set_xticklabels([])
            if r == 0:
                ax.set_title(CELL_NAME[cell], fontsize=9)
            if c == 0:
                ax.set_ylabel(f'{(rowlab or {}).get(r, col)}\nevents', fontsize=7)
            ax.text(0.98, 0.95, f'mixed mean {mixed.mean():.0f}\n({100 * mixed.mean() / max(hm, 1):.0f}% of human-mouse)',
                    transform=ax.transAxes, ha='right', va='top', fontsize=5.5)
    fig.suptitle(title, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(png, dpi=200)
    fig.savefig(png.replace('.png', '.pdf'))
    plt.close(fig)


def insert(x, sheet, png):
    wb = openpyxl.load_workbook(x)
    img = Image(png)
    img.width, img.height = img.width * 0.45, img.height * 0.45
    wb[sheet].add_image(img, 'K3')
    wb.save(x)


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else '5datasets'
    x = D + f'DSA_tables_{name}.xlsx'
    o2 = pd.read_csv(D + f'option2_splits_{name}.tsv', sep='\t')
    png2 = D + f'option2_splits_{name}.png'
    plot(o2, [(c, lambda b: b.index == b.index) for c in COLS],
         'Option 2: events separated by each split of the datasets (stage 1 FDR 5%); orange, human vs mouse; gray, '
         'the splits that mix the species; dashed line, mean of the mixed splits', png2)
    insert(x, 'option2_splits', png2)
    o1 = pd.read_csv(D + f'option1_splits_{name}.tsv', sep='\t')
    png1 = D + f'option1_splits_{name}.png'
    c1 = ['FDR 5%: significant', 'FDR 5%: significant + requirement', 'FDR 10%: significant',
          'FDR 10%: significant + requirement']
    plot(o1, [(c, lambda b: b.index == b.index) for c in c1],
         'Option 1: events significant in the Wilcoxon test of each split of the datasets (stage-1 events at the same FDR); '
         'orange, human vs mouse; gray, the splits that mix the species; dashed line, mean of the mixed splits', png1)
    insert(x, 'option1_splits', png1)
    print('plots added to', x)


if __name__ == '__main__':
    main()
