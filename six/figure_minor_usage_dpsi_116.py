#!/usr/bin/env python3
"""
Candidate figure: minor isoform usage vs the human-mouse difference of each event, per event type and immune cell type.
rows = event types (SE, A3SS, A5SS, AFE, ALE), columns = immune cell types;
x = minor isoform usage (%) of the event in the cell type (minor_usage_116.py: PSI of the less used isoform, averaged
over the six datasets); y = max|dPSI| of the event between species = the largest |mean of the human dataset PSIs -
mean of the mouse dataset PSIs| over the junctions of the event (dataset PSI = mean over the samples with >= 20 JSRs
in the event, as in the classification); skipped exons: the difference of the skipping junction. Testable events
only; colored by event type (complex events not shown).
Output: revision_6ds/figures/Figure_minor_usage_vs_dpsi.pdf / .tif and .tsv
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from figures_main_116 import CELL_LABEL, CELLS, ED, OUT, save

TYPES = {'SE': '#0072B2', 'A3SS': '#E69F00', 'A5SS': '#009E73', 'AFE': '#CC79A7', 'ALE': '#D55E00'}   # Okabe-Ito
NAMES = {'SE': 'Skipped exon', 'A3SS': "Alternative 3' splice site", 'A5SS': "Alternative 5' splice site",
         'AFE': 'Alternative first exon', 'ALE': 'Alternative last exon'}
MIN_READS = 20


def species_dpsi(counts, design, clu, cell):
    """max over the junctions of |mean human dataset PSI - mean mouse dataset PSI|, per event."""
    smp = design.index[design['cell'] == cell].intersection(counts.columns)
    tot = counts[smp].groupby(clu).transform('sum')
    psi = counts[smp] / tot.where(tot >= MIN_READS)
    means = pd.DataFrame({d: psi[design.index[(design['cell'] == cell) & (design['dataset'] == d)]
                                 .intersection(smp)].mean(axis=1) for d in ('H1', 'H2', 'H3', 'M1', 'M2', 'M3')})
    gap = (means[['H1', 'H2', 'H3']].mean(axis=1) - means[['M1', 'M2', 'M3']].mean(axis=1)).abs()
    return gap


def event_dpsi(gap, clu, se):
    """Per event: max over the junctions; skipped exons: the skipping junction (the longest intron), which carries
    the exclusion isoform alone, so that its PSI is on the same isoform scale as the minor isoform usage."""
    out = gap.groupby(clu).max()
    j = gap.index.to_series()
    length = j.str.split(':').str[2].astype(int) - j.str.split(':').str[1].astype(int)
    skip = length[clu.isin(se)].groupby(clu[clu.isin(se)]).idxmax()
    out.loc[skip.index] = gap.loc[skip.values].values
    return out


def main():
    counts = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', index_col=0)
    clu = pd.Series(counts.index.str.rsplit(':', n=1).str[1], index=counts.index)
    design = pd.read_csv(ED + 'sample_design.tsv', sep='\t').set_index('sample')
    mu = pd.read_csv(ED + 'revision/minor_usage_by_event.tsv', sep='\t',
                     usecols=['cell', 'cluster', 'event_type', 'status_high', 'minor_usage'])
    mu = mu[(mu['status_high'] != 'not testable') & mu['minor_usage'].notna() & mu['event_type'].isin(TYPES)]
    se = set(mu.loc[mu['event_type'] == 'SE', 'cluster'])
    gap = {cell: event_dpsi(species_dpsi(counts, design, clu, cell), clu, se) for cell in CELLS}
    mu['max_abs_dpsi_hm'] = [gap[c].get(k, np.nan) for c, k in zip(mu['cell'], mu['cluster'])]
    fig, axes = plt.subplots(len(TYPES), len(CELLS), figsize=(7.2, 7.6), sharex=True, sharey=True,
                             gridspec_kw={'wspace': 0.1, 'hspace': 0.18})
    for r, (t, c) in enumerate(TYPES.items()):
        for k, cell in enumerate(CELLS):
            ax = axes[r, k]
            x = mu[(mu['cell'] == cell) & (mu['event_type'] == t)]
            ax.scatter(100 * x['minor_usage'], x['max_abs_dpsi_hm'], s=3, color=c, alpha=0.5, lw=0, rasterized=True)
            ax.text(0.97, 0.95, f'n = {len(x):,}', transform=ax.transAxes, ha='right', va='top', fontsize=5.5)
            ax.set_xlim(-1, 51)
            ax.set_ylim(-0.02, 1.0)
            ax.tick_params(labelsize=5.5)
            if r == 0:
                ax.set_title(CELL_LABEL[cell], fontsize=7)
            if r == len(TYPES) - 1:
                ax.set_xlabel('minor isoform usage (%)', fontsize=6.5)
            if k == 0:
                ax.set_ylabel(f'{NAMES[t]}\nhuman-mouse ' + ('|ΔPSI| (skipping junction)' if t == 'SE' else 'max |ΔPSI|'),
                              fontsize=6)
    rows = [mu]
    save(fig, 'Figure_minor_usage_vs_dpsi')
    t = pd.concat(rows)
    t.to_csv(OUT + 'Figure_minor_usage_vs_dpsi.tsv', sep='\t', index=False)
    print(t.groupby('event_type')[['minor_usage', 'max_abs_dpsi_hm']].median().round(3).to_string())
    print('Spearman r:', round(t[['minor_usage', 'max_abs_dpsi_hm']].corr(method='spearman').iloc[0, 1], 3))


if __name__ == '__main__':
    main()
