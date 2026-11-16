#!/usr/bin/env python3
"""
Figure 1 of the revision. Panel A (orthologous exons / junctions schematic) is taken unchanged from
the submitted figure (generate_figure1_A4_v6.py output); panels B (datasets and comparisons) and
C (classification) are redrawn for the revised design, in the same style.
Usage: figure1_116.py <submitted figure 1 PDF>   -> revision_6ds/figures/Figure_1.pdf / .tif
"""
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations
from matplotlib.patches import FancyArrowPatch, Polygon
from PIL import Image

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'pdf.fonttype': 42})
OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_6ds/figures/'
BLUE, ORANGE, GREEN = '#2c5f8e', '#e8762c', '#6abf69'
INK = '#333333'
HUMAN_FC, HUMAN_EC, MOUSE_FC, MOUSE_EC = '#dbe7f3', '#2c5f8e', '#fbe3d2', '#b4541a'
W, H = 8.27, 11.69  # A4 inches


def tag(ax, x, y, text, species, fs=8):
    fc, ec = (HUMAN_FC, HUMAN_EC) if species == 'H' else (MOUSE_FC, MOUSE_EC)
    ax.text(x, y, text, fontsize=fs, ha='center', va='center', fontweight='bold', color=ec,
            bbox=dict(boxstyle='round,pad=0.25', fc=fc, ec=ec, lw=0.8))


def arrow(ax, p, q, rad=0.0, color=INK, lw=1.0, style='<|-|>'):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=7, color=color, lw=lw,
                                 connectionstyle=f'arc3,rad={rad}', shrinkA=2, shrinkB=2))


def diamond(ax, x, y, w, h, text, fs=9):
    ax.add_patch(Polygon([[x - w / 2, y], [x, y + h / 2], [x + w / 2, y], [x, y - h / 2]], closed=True,
                         fill=False, ec=INK, lw=1.2))
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, linespacing=1.5)


def down(ax, x0, y0, x1, y1, label=None, lx=0, ly=0):
    ax.plot([x0, x0, x1], [y0, y1 + 0.0, y1], color=INK, lw=1.2) if False else None
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=8, color=INK, lw=1.2,
                                 shrinkA=0, shrinkB=0))
    if label:
        ax.text(x0 + lx, (y0 + y1) / 2 + ly, label, fontsize=8.5, ha='center', va='center')


def elbow(ax, x0, y0, x1, y1, label=None, lpos=None):
    ax.plot([x0, x1], [y0, y0], color=INK, lw=1.2)
    ax.add_patch(FancyArrowPatch((x1, y0), (x1, y1), arrowstyle='-|>', mutation_scale=8, color=INK, lw=1.2,
                                 shrinkA=0, shrinkB=0))
    if label:
        ax.text(*(lpos or ((x0 + x1) / 2, y0 + 0.012)), label, fontsize=8.5, ha='center', va='bottom')


def psi_bars(fig, rect, profiles, title):
    """profiles: list of 5 (blue, green, orange) fractions for M2, M1, H1, H2, H3; 2 bars each."""
    ax = fig.add_axes(rect)
    x, ticks = 0, []
    for k, (lab, prof) in enumerate(zip(['M2', 'M1', 'H1', 'H2', 'H3'], profiles)):
        if k == 2:
            x += 0.6
        xs = [x, x + 0.9]
        ticks.append((np.mean(xs), lab))
        for xi, jitter in zip(xs, (0.0, 0.015)):
            b, g = prof[0] + jitter, prof[1] - jitter
            ax.bar(xi, b, color=BLUE, width=0.75)
            ax.bar(xi, g, bottom=b, color=GREEN, width=0.75)
            ax.bar(xi, 1 - b - g, bottom=b + g, color=ORANGE, width=0.75)
        if k not in (1, 4):
            ax.axvline(x + 1.6, color='#9e9e9e', lw=0.6, ls=(0, (2, 2)))
        x += 2.2
    ax.set_xticks([t for t, _ in ticks])
    ax.set_xticklabels([l for _, l in ticks], fontsize=8)
    ax.tick_params(axis='x', length=0)
    ax.set_yticks([])
    ax.set_ylim(0, 1)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.set_title(title, fontsize=8, pad=3)
    return ax


def dataset_pairs_panel(ax, label='B'):
    """Datasets and the dataset pairs compared one at a time (y in 0.74-0.99 of a 0.805 x A4 axes)."""
    ax.text(0.03, 0.985, label, fontsize=12, fontweight='bold', va='top')
    ax.text(0.03, 0.955, 'Immune\ncell types', fontsize=9, va='center')
    hx, mx = 0.30, 0.74
    hy = [0.955, 0.905, 0.855]
    my = [0.955, 0.905, 0.855]
    for y, (lab, acc) in zip(hy, [('Human H1', 'GSE60424'), ('Human H2', 'GSE115736'), ('Human H3', 'GSE118165')]):
        ax.text(hx, y + 0.008, lab, fontsize=9, fontweight='bold', ha='right', va='center')
        ax.text(hx, y - 0.012, acc, fontsize=8, ha='right', va='center')
    for y, lab in zip(my, ['Mouse M1', 'Mouse M2', 'Mouse M3']):
        ax.text(mx, y + 0.008, lab, fontsize=9, fontweight='bold', ha='left', va='center')
        ax.text(mx, y - 0.012, 'one study per cell type', fontsize=7, ha='left', va='center', color='#555555')
    for y1 in hy:
        for y2 in my:
            ax.plot([hx + 0.02, mx - 0.02], [y1, y2], color=ORANGE, lw=0.8, alpha=0.9)
    ax.text((hx + mx) / 2, 0.966, 'i', fontsize=9, style='italic', ha='center')
    for a, b, r in [(0, 1, 0.5), (1, 2, 0.5), (0, 2, 0.6)]:
        arrow(ax, (hx - 0.155, hy[a]), (hx - 0.155, hy[b]), rad=r, color='#2a78d6', lw=0.9)
        arrow(ax, (mx + 0.215, my[a]), (mx + 0.215, my[b]), rad=-r, color='#1baf7a', lw=0.9)
    ax.text(hx - 0.122, 0.905, 'ii', fontsize=9, style='italic', ha='center', va='center')
    ax.text(mx + 0.25, 0.905, 'iii', fontsize=9, style='italic', ha='center', va='center')
    ax.text(0.03, 0.79, 'Fibroblasts', fontsize=9, va='center')
    ax.text(hx, 0.798, 'Human', fontsize=9, fontweight='bold', ha='right', va='center')
    ax.text(hx, 0.778, 'E-MTAB-5919', fontsize=8, ha='right', va='center')
    ax.text(mx, 0.798, 'Mouse', fontsize=9, fontweight='bold', ha='left', va='center')
    ax.text(mx, 0.778, 'E-MTAB-5919', fontsize=8, ha='left', va='center')
    arrow(ax, (hx + 0.02, 0.79), (mx - 0.02, 0.79))
    ax.text((hx + mx) / 2, 0.772, 'iv', fontsize=9, style='italic', ha='center')
    ax.text(0.5, 0.745, 'i: human-mouse pairs (9; effect sizes only)   ii: human-human pairs (3)   iii: mouse-mouse pairs (3)   '
            'iv: fibroblasts', fontsize=7.5, ha='center', color='#555555')



def main():
    src = sys.argv[1]
    tmp = OUT + 'tmp_fig1'
    subprocess.run(['pdftoppm', '-r', '300', '-png', '-singlefile', src, tmp], check=True)
    page = Image.open(tmp + '.png').convert('RGB')
    panel_a = page.crop((0, 0, page.size[0], int(page.size[1] * 0.195)))  # panel A only (B starts at ~20%)

    fig = plt.figure(figsize=(W, H))
    axA = fig.add_axes([0, 1 - 0.195, 1, 0.195])
    axA.imshow(panel_a)
    axA.axis('off')
    ax = fig.add_axes([0, 0, 1, 1 - 0.195])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    # ---------------- panel B: human-mouse classification of each event (ANOVA -> Wilcoxon -> 2x rule)
    ax.text(0.03, 0.985, 'B', fontsize=12, fontweight='bold', va='top')
    ax.text(0.5, 0.975, 'Each immune cell type: three human and three mouse datasets (at least three samples each)',
            fontsize=9, ha='center', style='italic')
    for x, lab, acc in [(0.10, 'H1', 'GSE60424'), (0.20, 'H2', 'GSE115736'), (0.30, 'H3', 'GSE118165')]:
        tag(ax, x, 0.926, lab, 'H', fs=9.5)
        ax.text(x, 0.903, acc, fontsize=7, ha='center', va='top')
    ax.text(0.20, 0.950, 'Human (all five cell types)', fontsize=9, fontweight='bold', ha='center', color=HUMAN_EC)
    ax.text(0.70, 0.950, 'Mouse', fontsize=9, fontweight='bold', ha='center', color=MOUSE_EC)
    mouse = {'B': ('GSE287527', 'GSE307667', 'GSE231772'), 'T4': ('GSE215775', 'GSE184496', 'GSE134322'),
             'T8': ('GSE313371', 'GSE263132', 'GSE284441'), 'NK': ('GSE261448', 'GSE236819', 'GSE180053'),
             'Mo': ('GSE261448', 'GSE254202', 'PRJNA797236')}
    xs = [0.60, 0.72, 0.84]
    for x, lab in zip(xs, ['M1', 'M2', 'M3']):
        tag(ax, x, 0.926, lab, 'M', fs=9.5)
    for i, (cell, accs) in enumerate(mouse.items()):
        y = 0.902 - 0.0135 * i
        ax.text(0.525, y, cell, fontsize=7, ha='right', va='center', fontweight='bold')
        for x, acc in zip(xs, accs):
            ax.text(x, y, acc, fontsize=7, ha='center', va='center')
    ax.text(0.5, 0.815, 'ANOVA on the junction with the largest range of dataset PSIs (maximum minus minimum over '
            'the six datasets); range and gap in every junction', fontsize=8.5, ha='center', color=INK)
    ax.text(0.5, 0.795, 'PSI per sample (samples with ≥ 20 JSRs in the event); dataset PSI = mean over its samples',
            fontsize=7.5, ha='center', color='#555555')
    diamond(ax, 0.5, 0.74, 0.44, 0.065, 'PSI in at least two samples\nof each of the six datasets?', fs=8.5)
    elbow(ax, 0.5 - 0.22, 0.74, 0.10, 0.70, 'No', (0.24, 0.743))
    ax.text(0.10, 0.693, 'Not testable', fontsize=8.5, ha='center', va='top')
    down(ax, 0.5, 0.7075, 0.5, 0.67, 'Yes', lx=0.03)
    diamond(ax, 0.5, 0.63, 0.62, 0.08, 'ANOVA over the six datasets\nP > 0.1 and dataset PSIs\nwithin a range of 0.05?', fs=8)
    elbow(ax, 0.5 + 0.31, 0.63, 0.88, 0.59, 'Yes', (0.84, 0.633))
    ax.text(0.88, 0.583, 'Splicing\nconserved', fontsize=9.5, fontweight='bold', ha='center', va='top', color='#2171b5')
    down(ax, 0.5, 0.59, 0.5, 0.55, 'No', lx=0.03)
    diamond(ax, 0.5, 0.49, 0.70, 0.12, 'ANOVA adjusted P ≤ 0.05, and in a junction:\n|mean human - mean mouse PSI| ≥ 0.1\nand smallest human-mouse |ΔPSI| >\n1.5 × largest within-species |ΔPSI|?', fs=8)
    elbow(ax, 0.5 - 0.35, 0.49, 0.10, 0.425, 'Yes', (0.13, 0.493))
    ax.text(0.10, 0.418, 'Differentially\nspliced', fontsize=9.5, fontweight='bold', ha='center', va='top',
            color='#238b45')
    elbow(ax, 0.5 + 0.35, 0.49, 0.90, 0.425, 'No', (0.87, 0.493))
    ax.text(0.90, 0.418, 'Inconclusive', fontsize=9.5, fontweight='bold', ha='center', va='top', color='#757575')
    ax.text(0.5, 0.37, '|ΔPSI| of the junction between dataset PSIs: three human-human, three mouse-mouse and nine '
            'human-mouse dataset pairs', fontsize=7.5, ha='center', color='#555555')
    ax.text(0.5, 0.348, 'Single dataset pairs (human-human, mouse-mouse), fibroblasts and sex, age and tissue '
            'comparisons: Figure S1', fontsize=7.5, ha='center', color='#555555')
    ax.text(0.03, 0.03, 'Figure 1', fontsize=10, fontweight='bold')
    fig.savefig(OUT + 'Figure_1.pdf')
    fig.savefig(OUT + 'Figure_1.tif', dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
    print('wrote Figure_1')


if __name__ == '__main__':
    main()
