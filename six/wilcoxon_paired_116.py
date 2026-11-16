#!/usr/bin/env python3
"""
Paired Wilcoxon signed-rank tests of effect sizes on the same events (not part of the manuscript; separate folder).

Effect size of an event in a comparison = max|dPSI| over its junctions (LeafCutter fitted PSI, effects/event_effects).
Human-mouse effect of an event = mean max|dPSI| over the six human-mouse dataset pairs of the cell type.

Sheet within_species: per cell type, the events testable in all ten dataset pairs; human-mouse effect vs each
  within-species pair (H1-H2, H1-H3, H2-H3, M1-M2) and vs the mean of the four.
Sheet pair_by_pair: every human-mouse pair vs every within-species pair, on the same events.
Sheet controls: per sex/age/tissue control comparison, the events testable in the control and in all six
  human-mouse pairs of its cell type; human-mouse effect vs the control.
Test: one-sided paired Wilcoxon signed-rank test, alternative = human-mouse effect larger. Also reported: medians,
the median of the per-event differences and the % of events with a larger human-mouse effect.

Output: ensembl116/events/id50_wilcoxon/paired_wilcoxon.(tsv per sheet|xlsx)
"""
import itertools
import os

import numpy as np
import pandas as pd
from scipy import stats

import build_supp_tables_116 as bst

E = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50_6ds/'
OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50_wilcoxon/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
HM = [f'{h}_{m}' for h in ('H1', 'H2', 'H3') for m in ('M1', 'M2')]
WITHIN = ['H1_H2', 'H1_H3', 'H2_H3', 'M1_M2']
CONTROL_LABEL = {'human_sex': 'Human sex', 'human_age': 'Human age', 'human_tissue': 'Human tissue',
                 'mouse_sex': 'Mouse sex', 'mouse_age': 'Mouse age', 'mouse_tissue': 'Mouse tissue (NK, spleen vs liver)',
                 'mouse_tissue_m1': 'Mouse tissue (Mo, bone marrow vs blood)',
                 'mouse_tissue_mono': 'Mouse tissue (Mo, blood vs spleen)'}


def wide(path):
    e = pd.read_csv(path, sep='\t', usecols=['cluster', 'comparison', 'cell', 'testable', 'max_abs_dpsi'])
    e = e[e['testable'] & e['max_abs_dpsi'].notna()]
    return e.pivot_table(index=['cell', 'cluster'], columns='comparison', values='max_abs_dpsi')


def test(a, b):
    d = a - b
    return {'events': len(d), 'median_human_mouse': round(np.median(a), 4), 'median_other': round(np.median(b), 4),
            'median_difference': round(np.median(d), 4), 'pct_events_human_mouse_larger': round(100 * (d > 0).mean(), 1),
            'wilcoxon_p': stats.wilcoxon(a, b, alternative='greater').pvalue}


def main():
    os.makedirs(OUT, exist_ok=True)
    w = wide(E + 'effects/event_effects.tsv')
    within, pairs, controls = [], [], []
    for cell in CELLS:
        d = w.loc[cell, HM + WITHIN].dropna()
        hm = d[HM].mean(axis=1)
        for o in WITHIN + ['mean of the four']:
            other = d[WITHIN].mean(axis=1) if o == 'mean of the four' else d[o]
            within.append({'cell': cell, 'compared_with': o.replace('_', '-'), **test(hm, other)})
        for h, o in itertools.product(HM, WITHIN):
            pairs.append({'cell': cell, 'human_mouse_pair': h.replace('_', '-'), 'within_pair': o.replace('_', '-'),
                          **test(d[h], d[o])})
    c = wide(E + 'effects_confounders/event_effects.tsv')
    for (cell, comp) in sorted({(i[0], col) for i in c.index for col in c.columns}):
        if comp not in CONTROL_LABEL or cell not in CELLS:
            continue
        x = pd.concat([w.loc[cell, HM], c.loc[cell, [comp]]], axis=1, join='inner').dropna()
        if len(x) == 0:
            continue
        controls.append({'control': CONTROL_LABEL[comp], 'cell': cell, **test(x[HM].mean(axis=1), x[comp])})
    sheets = {'within_species': pd.DataFrame(within), 'pair_by_pair': pd.DataFrame(pairs),
              'controls': pd.DataFrame(controls)}
    for k, v in sheets.items():
        v.to_csv(OUT + f'paired_wilcoxon_{k}.tsv', sep='\t', index=False)
    bst.OUT = OUT
    bst.write('paired_wilcoxon', 'Paired Wilcoxon tests: human-mouse vs within-species and control effect sizes',
              'Effect size = max|dPSI| of the event in a comparison; human-mouse = mean over the six human-mouse dataset '
              'pairs of the cell type. Same events in each test (testable in all compared comparisons). One-sided '
              'paired Wilcoxon signed-rank test, alternative: human-mouse effect larger.', sheets)
    pd.set_option('display.width', 200)
    for k, v in sheets.items():
        print('==', k)
        print(v.head(40).to_string(index=False))


if __name__ == '__main__':
    main()
