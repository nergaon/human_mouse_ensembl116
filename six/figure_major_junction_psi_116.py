#!/usr/bin/env python3
"""
Candidate figure: PSI of the human major junction of each event, in human and in mouse, per cell type.

For each event (LeafCutter cluster) and cell type, the human major junction is the junction with the highest human
PSI (mean of the dataset PSIs of H1, H2 and H3; fibroblasts: FH). Its PSI is shown in human (same mean) and in
mouse (mean of M1 and M2; fibroblasts: FM). Dataset PSI = mean over the samples with >= 20 reads in the event (as in
the classification). Events: immune cell types, the events with a PSI in all five datasets (testable in the
classification); fibroblasts, the events with a PSI in FH and FM.

Output: ensembl116/revision_6ds/candidate_figures/Figure_major_junction_psi.(pdf|png), major_junction_psi.tsv
"""
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import classify_separation_116 as cs  # noqa: E402

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})
E = cs.BASE + 'events/id50_6ds/'
OUT = cs.BASE + 'revision_6ds/candidate_figures/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono', 'Fibro']
CELL_NAME = {'CD4T': 'CD4 T', 'CD8T': 'CD8 T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Monocytes', 'Fibro': 'Fibroblasts'}
SPECIES_COLOR = {'Human': '#2c5f8e', 'Mouse': '#b4541a'}


def dataset_psi(counts, design, cell, datasets):
    cluster = counts.index.str.split(':').str[3]
    psi = {}
    for ds in datasets:
        s = design[(design['cell'] == cell) & (design['dataset'] == ds)].index.intersection(counts.columns)
        tot = counts[s].groupby(cluster).transform('sum')
        psi[ds] = (counts[s] / tot.where(tot >= cs.MIN_READS)).mean(axis=1)
    psi = pd.DataFrame(psi)
    ok = psi.notna().all(axis=1).groupby(cluster).transform('all')
    psi = psi[ok.to_numpy()]
    return psi, psi.index.str.split(':').str[3]


def main():
    counts = pd.read_csv(E + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    rows = []
    for cell in CELLS:
        hum, mou = (['FH'], ['FM']) if cell == 'Fibro' else (['H1', 'H2', 'H3'], ['M1', 'M2'])
        psi, clu = dataset_psi(counts, design, cell, hum + mou)
        d = pd.DataFrame({'cluster': clu, 'junction': psi.index, 'human': psi[hum].mean(axis=1).to_numpy(),
                          'mouse': psi[mou].mean(axis=1).to_numpy()})
        major = d.loc[d.groupby('cluster')['human'].idxmax()]
        rows.append(major.assign(cell=cell))
    t = pd.concat(rows, ignore_index=True)
    t.to_csv(OUT + 'major_junction_psi.tsv', sep='\t', index=False, float_format='%.4f')

    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    pos, labels, k = [], [], 0
    for cell in CELLS:
        d = t[t['cell'] == cell]
        for sp, col in (('human', SPECIES_COLOR['Human']), ('mouse', SPECIES_COLOR['Mouse'])):
            v = d[sp].to_numpy()
            parts = ax.violinplot(v, positions=[k], widths=0.8, showextrema=False, bw_method=0.15)
            for b in parts['bodies']:
                b.set_facecolor(col)
                b.set_edgecolor('none')
                b.set_alpha(0.45)
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.plot([k, k], [q1, q3], color=col, lw=3, solid_capstyle='butt')
            ax.plot(k, med, 'o', color='white', ms=3, mec=col, mew=0.8)
            pos.append(k)
            labels.append('Human' if sp == 'human' else 'Mouse')
            k += 1
        ax.text(k - 1.5, 1.07, f'n = {len(d):,}', ha='center', fontsize=6, color='#555555')
        ax.text(k - 1.5, -0.17, CELL_NAME[cell], ha='center', va='top', fontsize=7.5, fontweight='bold',
                transform=ax.get_xaxis_transform())
        k += 0.6
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.set_ylim(0, 1.12)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_ylabel('PSI of the human major junction')
    ax.set_title('PSI of the junction with the highest human PSI in each event, in human and in mouse', fontsize=8,
                 loc='left')
    fig.text(0.01, 0.005, 'Human: mean PSI of H1, H2, H3 (fibroblasts: FH); mouse: mean of M1, M2 (fibroblasts: FM). '
             'Violin: distribution over events; bar: interquartile range; circle: median. n: events.', fontsize=6)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    for ext in ('pdf', 'png'):
        fig.savefig(OUT + f'Figure_major_junction_psi.{ext}', dpi=300)
    plt.close(fig)
    s = t.groupby('cell')[['human', 'mouse']].describe(percentiles=[0.1, 0.25, 0.5]).round(2)
    print(s.loc[CELLS, [(a, b) for a in ('human', 'mouse') for b in ('count', '10%', '25%', '50%')]])


if __name__ == '__main__':
    main()
