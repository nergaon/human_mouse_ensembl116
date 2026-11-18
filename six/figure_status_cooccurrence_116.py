#!/usr/bin/env python3
"""
Candidate figure: co-occurrence of the human-mouse classes of an event across the immune cell types.
Rows and columns = (class, cell type) for differentially spliced (DS), splicing conserved (SC), inconclusive and not
testable x the five immune cell types; value = number of events with the row class in the row cell type and the
column class in the column cell type (diagonal of each block: number of events with that class in that cell type).
Color capped at 200 events so that the not-testable blocks do not dominate; every cell shows its number.
Output: revision_6ds/figures/Figure_status_cooccurrence.pdf / .tif and .tsv
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figures_main_116 import CELL_LABEL, CELLS, ED, OUT, save

CLASSES = [('divergent', 'DS'), ('conserved', 'SC'), ('inconclusive', 'Inc'), ('not testable', 'NT')]
CAP = 200


def main():
    st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t', usecols=['cell', 'cluster', 'status_high'])
    w = st.pivot(index='cluster', columns='cell', values='status_high')[CELLS]
    keys = [(c, cell) for c, _ in CLASSES for cell in CELLS]
    ind = pd.DataFrame({k: (w[k[1]] == k[0]).astype(int) for k in keys})
    m = ind.T.dot(ind)                                    # events with both (class, cell) pairs
    labels = [f'{short} {CELL_LABEL[cell]}' for (c, cell), short in
              zip(keys, [s for _, s in CLASSES for _ in CELLS])]
    m.index = m.columns = labels
    m.to_csv(OUT + 'Figure_status_cooccurrence.tsv', sep='\t')
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    ax.imshow(np.minimum(m.values, CAP), cmap='Greys', vmin=0, vmax=CAP, interpolation='none')
    for i in range(len(m)):
        for j in range(len(m)):
            v = m.values[i, j]
            ax.text(j, i, f'{v:,}' if v < 10000 else f'{v / 1000:.0f}k', ha='center', va='center', fontsize=4.2,
                    color='white' if min(v, CAP) > CAP * 0.55 else '#1f2328')
    for b in range(1, len(CLASSES)):                      # block borders between classes
        ax.axhline(b * len(CELLS) - 0.5, color='#c0392b', lw=0.8)
        ax.axvline(b * len(CELLS) - 0.5, color='#c0392b', lw=0.8)
    ax.set_xticks(range(len(m)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_yticks(range(len(m)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    sm = matplotlib.cm.ScalarMappable(cmap='Greys', norm=matplotlib.colors.Normalize(0, CAP))
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label(f'events (color capped at {CAP})', fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    save(fig, 'Figure_status_cooccurrence')
    ds_sc = m.loc[[l for l in labels if l.startswith('DS')], [l for l in labels if l.startswith('SC')]]
    print('DS x SC (events DS in the row cell type and SC in the column cell type):')
    print(ds_sc.to_string())


if __name__ == '__main__':
    main()
