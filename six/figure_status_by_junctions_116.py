#!/usr/bin/env python3
"""
Candidate figure: human-mouse status of the events by their number of junctions. One panel per immune cell type;
x = number of junctions of the event (2, 3, 4, 5, 6, >= 7), each bar = 100% of the events with that number of
junctions (their number above the bar), stacked by status (differentially spliced, splicing conserved, inconclusive, not testable; Figure 1B).
Output: revision_6ds/figures/Figure_status_by_junctions.pdf / .tif and the numbers as .tsv
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from figures_main_116 import CELL_LABEL, CELLS, ED, HC_COLORS, NOT_TESTABLE, OUT, save

ORDER = ['divergent', 'conserved', 'inconclusive', 'not testable']          # bottom to top
COLORS = {**HC_COLORS, 'not testable': NOT_TESTABLE}
NAMES = {'divergent': 'Differentially spliced', 'conserved': 'Splicing conserved', 'inconclusive': 'Inconclusive',
         'not testable': 'Not testable'}
MAX_BIN = 7


def main():
    st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t', usecols=['cell', 'cluster', 'status_high'])
    ids = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', usecols=[0]).iloc[:, 0]
    nj = ids.str.rsplit(':', n=1).str[1].value_counts()
    st['junctions'] = st['cluster'].map(nj).clip(upper=MAX_BIN)
    assert st['junctions'].notna().all()
    bins = list(range(2, MAX_BIN + 1))
    ticks = [str(b) for b in bins[:-1]] + [f'≥{MAX_BIN}']
    n = st.groupby(['cell', 'junctions']).size()
    pct = st.groupby(['cell', 'junctions'])['status_high'].value_counts(normalize=True).mul(100).unstack(fill_value=0)
    pct = pct.reindex(columns=ORDER, fill_value=0)
    out = pct.round(1).join(n.rename('events')).reset_index()
    out.to_csv(OUT + 'Figure_status_by_junctions.tsv', sep='\t', index=False)

    fig, axes = plt.subplots(1, len(CELLS), figsize=(7.2, 2.6), sharey=True, gridspec_kw={'wspace': 0.12})
    x = np.arange(len(bins))
    for i, (ax, cell) in enumerate(zip(axes, CELLS)):
        bottom = np.zeros(len(bins))
        for s in ORDER:
            v = pct.loc[cell].reindex(bins)[s].to_numpy()
            ax.bar(x, v, bottom=bottom, width=0.75, color=COLORS[s], edgecolor='#9a9a9a' if s == 'not testable'
                   else 'white', linewidth=0.4)
            bottom += v
        for xi, b in zip(x, bins):
            ax.text(xi, 101, f'{n.loc[(cell, b)]:,}', ha='center', va='bottom', fontsize=4.5, rotation=90)
        ax.set_xticks(x, ticks)
        ax.set_xlabel('Junctions in the event')
        ax.set_ylim(0, 100)
        ax.set_title(CELL_LABEL[cell], fontsize=7, pad=16)
        if i == 0:
            ax.set_ylabel('% of events')
            ax.text(-0.7, 101, 'events:', ha='right', va='bottom', fontsize=4.5)
    fig.legend([Patch(facecolor=COLORS[s], edgecolor='#9a9a9a' if s == 'not testable' else 'none') for s in ORDER[::-1]],
               [NAMES[s] for s in ORDER[::-1]], loc='center left', bbox_to_anchor=(0.91, 0.5), frameon=False,
               fontsize=6)
    save(fig, 'Figure_status_by_junctions')
    print(out[out['cell'] == 'NK'].to_string(index=False))


if __name__ == '__main__':
    main()
