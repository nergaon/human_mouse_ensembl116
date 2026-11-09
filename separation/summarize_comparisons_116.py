#!/usr/bin/env python3
"""
Summary of the effect-size / equivalence classification per comparison type and cell type.

Comparison types: cross (human-mouse dataset pairs), HH, MM (within species), fibro (FH-FM);
'matched_' prefix = M1 restricted to naive T cells and blood monocytes.
Input:  <effects_dir>/event_effects.tsv, psi_correlations.tsv (effect_size_analysis.py)
Output: <effects_dir>/comparison_summary.tsv (per comparison) and comparison_type_summary.tsv
Usage:  summarize_comparisons_116.py <effects_dir>
"""
import os
import sys

import pandas as pd

MARGINS = ('0.025', '0.05', '0.1')


def kind(comparison):
    if comparison.startswith('matched_'):
        return 'matched_' + kind(comparison[len('matched_'):])
    if comparison == 'FH_FM':
        return 'fibro'
    a, b = comparison.split('_')
    if a[0] != b[0]:
        return 'cross'
    return 'HH' if a[0] == 'H' else 'MM'


def summarize(events, corr):
    rows = []
    for (comp, cell), g in events.groupby(['comparison', 'cell']):
        g = g[g['status_m0.1'] != 'not_testable']
        row = dict(comparison=comp, kind=kind(comp), cell=cell, testable=len(g),
                   pct_DS=100 * (g['status_m0.1'] == 'DS').mean(),
                   median_max_abs_dpsi=g['max_abs_dpsi'].median())
        for m in MARGINS:
            row[f'pct_equivalent_m{m}'] = 100 * (g[f'status_m{m}'] == 'equivalent').mean()
        rows.append(row)
    per_comp = pd.DataFrame(rows).merge(corr[['comparison', 'cell', 'pearson']], on=['comparison', 'cell'])
    per_kind = per_comp.groupby(['kind', 'cell']).agg(
        comparisons=('comparison', 'nunique'), testable_mean=('testable', 'mean'),
        pct_DS_min=('pct_DS', 'min'), pct_DS_max=('pct_DS', 'max'),
        pearson_min=('pearson', 'min'), pearson_max=('pearson', 'max'),
        median_abs_dpsi=('median_max_abs_dpsi', 'median'),
        **{f'pct_eq_m{m}_mean': (f'pct_equivalent_m{m}', 'mean') for m in MARGINS}).reset_index()
    return per_comp, per_kind


def main():
    d = sys.argv[1]
    events = pd.read_csv(os.path.join(d, 'event_effects.tsv'), sep='\t')
    corr = pd.read_csv(os.path.join(d, 'psi_correlations.tsv'), sep='\t')
    per_comp, per_kind = summarize(events, corr)
    per_comp.round(3).to_csv(os.path.join(d, 'comparison_summary.tsv'), sep='\t', index=False)
    per_kind.round(3).to_csv(os.path.join(d, 'comparison_type_summary.tsv'), sep='\t', index=False)
    print(per_kind.round(2).to_string(index=False))


if __name__ == '__main__':
    main()
