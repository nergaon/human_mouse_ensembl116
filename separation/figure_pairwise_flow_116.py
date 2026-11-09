#!/usr/bin/env python3
"""
Figure S1 of the revision: the status of an alternative splicing event in a single comparison of two datasets
(LeafCutter p.adjust, max|ΔPSI| and the 90% CI of ΔPSI), used for the within-species, single human-mouse pair,
fibroblast and sex / age / tissue comparisons. The human-mouse classification across all five datasets is in
Figure 1C.  Usage: figure_s11_116.py   -> revision_separation/figures/Figure_S1.pdf / .tif
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from figure1_116 import OUT, INK, H as A4_H, dataset_pairs_panel, diamond, down, elbow

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'pdf.fonttype': 42})
W, H_LOW = 8.27, 6.6
H_TOP = 0.25 * 0.805 * A4_H   # panel A at the same scale as in the submitted Figure 1B
H = H_LOW + H_TOP + 0.1

COMPARISONS = [
    ('Human-human', '3 dataset pairs (H1-H2, H1-H3, H2-H3), five immune cell types', 'Figure 2A, Table 1'),
    ('Mouse-mouse', '1 dataset pair (M1-M2), five immune cell types', 'Figure 2B, Table 1'),
    ('Fibroblasts', 'human vs mouse dermal fibroblasts (E-MTAB-5919)', 'Table 1, Figure 6'),
    ('Sex', 'male vs female (human, five cell types; mouse B cells)', 'Figure S7, Table S12'),
    ('Age', 'young vs old (human and mouse naive CD4 T cells)', 'Figure S7, Table S12'),
    ('Tissue', 'human NK blood vs spleen; mouse NK spleen vs liver;\nmouse monocytes bone marrow vs blood',
     'Figure S7, Table S12'),
]


def main():
    fig = plt.figure(figsize=(W, H))
    top = fig.add_axes([0, (H - H_TOP - 0.05) / H, 1, H_TOP / H])
    top.set_xlim(0, 1)
    top.set_ylim(0.735, 0.995)
    top.axis('off')
    dataset_pairs_panel(top, label='A')
    ax = fig.add_axes([0, 0, 1, H_LOW / H])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    # ---------------- panel A: the comparisons that use the single-comparison status
    ax.text(0.03, 0.975, 'B', fontsize=12, fontweight='bold', va='top')
    ax.text(0.5, 0.955, 'Comparisons of two datasets (or two sample groups) in one cell type', fontsize=9,
            ha='center', style='italic')
    y = 0.91
    for name, what, where in COMPARISONS:
        ax.text(0.07, y, name, fontsize=8.5, fontweight='bold', va='center')
        ax.text(0.30, y, what, fontsize=8.5, va='center', linespacing=1.3)
        ax.text(0.95, y, where, fontsize=8, va='center', ha='right', color='#555555')
        y -= 0.05 if '\n' not in what else 0.065
    ax.text(0.5, y + 0.005, 'The human-mouse classification of each event across all five datasets is shown in '
            'Figure 1B, C.', fontsize=7.5, ha='center', color='#555555')
    ax.plot([0.06, 0.94], [y - 0.02, y - 0.02], color='#bdbdbd', lw=0.8)
    # ---------------- panel B: decision process (same layout as the submitted Figure 1C)
    t = y - 0.045
    ax.text(0.03, t + 0.01, 'C', fontsize=12, fontweight='bold', va='top')
    d1 = t - 0.06
    diamond(ax, 0.5, d1, 0.46, 0.085, 'Testable?\n(≥ g samples with ≥ 20 JSRs\nin each group)', fs=8)
    elbow(ax, 0.5 - 0.23, d1, 0.10, d1 - 0.055, 'No', (0.22, d1 + 0.004))
    ax.text(0.10, d1 - 0.065, 'Not testable', fontsize=8.5, ha='center', va='top')
    d2 = d1 - 0.12
    down(ax, 0.5, d1 - 0.0425, 0.5, d2 + 0.0375, 'Yes', lx=0.03)
    diamond(ax, 0.5, d2, 0.36, 0.075, 'p.adjust ≤ 0.05 and\nmax|ΔPSI| ≥ 0.1?', fs=8.5)
    elbow(ax, 0.5 - 0.18, d2, 0.20, d2 - 0.06, 'Yes', (0.26, d2 + 0.004))
    ax.text(0.20, d2 - 0.07, 'Differentially\nspliced (DS)', fontsize=8.5, ha='center', va='top', color='#238b45')
    d3 = d2 - 0.085
    elbow(ax, 0.5 + 0.18, d2, 0.80, d3 + 0.0375, 'No', (0.74, d2 + 0.004))
    diamond(ax, 0.80, d3, 0.30, 0.075, '90% CI of ΔPSI of every\njunction within ±0.1?', fs=8)
    elbow(ax, 0.80 - 0.15, d3, 0.58, d3 - 0.06, 'Yes', (0.62, d3 + 0.004))
    ax.text(0.58, d3 - 0.07, 'Equivalent', fontsize=8.5, ha='center', va='top', color='#2171b5')
    down(ax, 0.80, d3 - 0.0375, 0.80, d3 - 0.06, None)
    ax.text(0.83, d3 - 0.048, 'No', fontsize=8.5, ha='left', va='center')
    ax.text(0.80, d3 - 0.07, 'Inconclusive', fontsize=8.5, ha='center', va='top', color='#757575')
    ax.text(0.03, 0.02, 'Figure S1', fontsize=10, fontweight='bold')
    fig.savefig(OUT + 'Figure_S1.pdf')
    fig.savefig(OUT + 'Figure_S1.tif', dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
    print('wrote Figure_S1')


if __name__ == '__main__':
    main()
