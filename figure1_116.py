#!/usr/bin/env python3
"""
Figure 1 of the revision. Panel A (orthologous exons / junctions schematic) is taken unchanged from
the submitted figure (generate_figure1_A4_v6.py output); panels B (datasets and comparisons) and
C (classification) are redrawn for the revised design, in the same style.
Usage: figure1_116.py <submitted figure 1 PDF>   -> revision_files/figures/Figure_1.pdf / .tif
"""
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon
from PIL import Image

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'pdf.fonttype': 42})
OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_files/figures/'
BLUE, ORANGE, GREEN = '#2c5f8e', '#e8762c', '#6abf69'
INK = '#333333'
W, H = 8.27, 11.69  # A4 inches


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
    # ---------------- panel B: datasets and comparisons
    ax.text(0.03, 0.985, 'B', fontsize=12, fontweight='bold', va='top')
    ax.text(0.03, 0.955, 'Immune\ncell types', fontsize=9, va='center')
    hx, mx = 0.30, 0.74
    hy = [0.955, 0.905, 0.855]
    my = [0.93, 0.88]
    for y, (lab, acc) in zip(hy, [('Human H1', 'GSE60424'), ('Human H2', 'GSE115736'), ('Human H3', 'GSE118165')]):
        ax.text(hx, y + 0.008, lab, fontsize=9, fontweight='bold', ha='right', va='center')
        ax.text(hx, y - 0.012, acc, fontsize=8, ha='right', va='center')
    for y, (lab, acc) in zip(my, [('Mouse M1', 'GSE116177'),
                                  ('Mouse M2', 'GSE287527, GSE261448,\nENCODE4, GSE313371')]):
        ax.text(mx, y + 0.008, lab, fontsize=9, fontweight='bold', ha='left', va='center')
        ax.text(mx, y - 0.016 if '\n' in acc else y - 0.012, acc, fontsize=8, ha='left', va='center',
                linespacing=1.2)
    for y1 in hy:
        for y2 in my:
            ax.plot([hx + 0.02, mx - 0.02], [y1, y2], color=ORANGE, lw=0.8, alpha=0.9)
    ax.text((hx + mx) / 2, 0.966, 'i', fontsize=9, style='italic', ha='center')
    for a, b, r in [(0, 1, 0.5), (1, 2, 0.5), (0, 2, 0.6)]:
        arrow(ax, (hx - 0.155, hy[a]), (hx - 0.155, hy[b]), rad=r, color='#2a78d6', lw=0.9)
    ax.text(hx - 0.122, 0.905, 'ii', fontsize=9, style='italic', ha='center', va='center')
    arrow(ax, (mx + 0.215, my[0] + 0.004), (mx + 0.215, my[1] + 0.006), rad=-0.6, color='#1baf7a', lw=0.9)
    ax.text(mx + 0.25, 0.905, 'iii', fontsize=9, style='italic', ha='center', va='center')
    ax.text(0.03, 0.79, 'Fibroblasts', fontsize=9, va='center')
    ax.text(hx, 0.798, 'Human', fontsize=9, fontweight='bold', ha='right', va='center')
    ax.text(hx, 0.778, 'E-MTAB-5919', fontsize=8, ha='right', va='center')
    ax.text(mx, 0.798, 'Mouse', fontsize=9, fontweight='bold', ha='left', va='center')
    ax.text(mx, 0.778, 'E-MTAB-5919', fontsize=8, ha='left', va='center')
    arrow(ax, (hx + 0.02, 0.79), (mx - 0.02, 0.79))
    ax.text((hx + mx) / 2, 0.772, 'iv', fontsize=9, style='italic', ha='center')
    ax.text(0.5, 0.745, 'i: human-mouse pairs (6)   ii: human-human pairs (3)   iii: mouse-mouse pair (1)   '
            'iv: fibroblasts', fontsize=7.5, ha='center', color='#555555')

    # ---------------- panel C: classification
    ax.text(0.03, 0.715, 'C', fontsize=12, fontweight='bold', va='top')
    ax.text(0.5, 0.70, 'Each comparison (dataset pair, cell type)', fontsize=8.5, ha='center', style='italic')
    diamond(ax, 0.5, 0.655, 0.27, 0.06, 'Comparable?')
    elbow(ax, 0.5 - 0.135, 0.655, 0.10, 0.615, 'No', (0.25, 0.658))
    ax.text(0.10, 0.60, 'Not testable', fontsize=8.5, ha='center', va='top')
    down(ax, 0.5, 0.625, 0.5, 0.592, 'Yes', lx=0.03)
    diamond(ax, 0.5, 0.562, 0.36, 0.06, 'p.adjust ≤ 0.05 and\nmax|ΔPSI| ≥ 0.1?', fs=8.5)
    elbow(ax, 0.5 - 0.18, 0.562, 0.20, 0.515, 'Yes', (0.26, 0.565))
    ax.text(0.20, 0.51, 'Differentially\nspliced (DS)', fontsize=8.5, ha='center', va='top', color='#238b45')
    elbow(ax, 0.5 + 0.18, 0.562, 0.80, 0.53, 'No', (0.74, 0.565))
    diamond(ax, 0.80, 0.50, 0.30, 0.06, '90% CI of ΔPSI\nwithin ±0.1?', fs=8.5)
    elbow(ax, 0.80 - 0.15, 0.50, 0.58, 0.455, 'Yes', (0.62, 0.503))
    ax.text(0.58, 0.45, 'Equivalent', fontsize=8.5, ha='center', va='top', color='#2171b5')
    down(ax, 0.80, 0.47, 0.80, 0.455, None)
    ax.text(0.83, 0.462, 'No', fontsize=8.5, ha='left', va='center')
    ax.text(0.80, 0.45, 'Inconclusive', fontsize=8.5, ha='center', va='top', color='#757575')
    # high confidence
    ax.plot([0.06, 0.94], [0.405, 0.405], color='#bdbdbd', lw=0.8)
    ax.text(0.5, 0.39, 'Each event and cell type, across all comparisons', fontsize=8.5, ha='center',
            style='italic')
    boxes = [(0.19, 'Differentially spliced', '#238b45',
              'DS in two human-mouse pairs\nwith no shared dataset (i);\nequivalent in none (i)'),
             (0.5, 'Equivalent', '#2171b5', 'Equivalent in two human-mouse\npairs with no shared dataset (i);\nDS in none (i)'),
             (0.81, 'Inconclusive', '#757575', 'All other\ntestable events')]
    for x, title, color, text in boxes:
        ax.text(x, 0.365, title, fontsize=9.5, fontweight='bold', ha='center', va='top', color=color)
        ax.text(x, 0.338, text, fontsize=8, ha='center', va='top', linespacing=1.35)
    ax.text(0.5, 0.274, 'Confidence of differentially spliced and equivalent events:  high – no DS in ii, iii, and at least one '
            'human-human pair (ii) and the mouse-mouse pair (iii) equivalent;  low – all others',
            fontsize=6.3, ha='center', va='top', color='#333333')
    psi_bars(fig, [0.07, 0.06, 0.25, 0.11],
             [(0.10, 0.65), (0.11, 0.64), (0.45, 0.25), (0.44, 0.27), (0.46, 0.24)], 'Mouse ≠ Human')
    psi_bars(fig, [0.38, 0.06, 0.25, 0.11],
             [(0.30, 0.42), (0.31, 0.41), (0.30, 0.43), (0.29, 0.42), (0.31, 0.42)], 'Mouse = Human')
    psi_bars(fig, [0.69, 0.06, 0.25, 0.11],
             [(0.10, 0.70), (0.40, 0.35), (0.20, 0.55), (0.45, 0.20), (0.30, 0.45)], 'Datasets disagree')
    ax.text(0.03, 0.03, 'Figure 1', fontsize=10, fontweight='bold')
    fig.savefig(OUT + 'Figure_1.pdf')
    fig.savefig(OUT + 'Figure_1.tif', dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
    print('wrote Figure_1')


if __name__ == '__main__':
    main()
