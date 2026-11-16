#!/usr/bin/env python3
"""
Single human-mouse dataset pairs (not part of the manuscript; separate folder): for each of the six human-mouse pairs
and each immune cell type, the number of events tested by LeafCutter, the number significant (p.adjust <= 0.05) and
significant with max|dPSI| >= 0.1, and their FDR.

Significant: FDR (LeafCutter p.adjust, Benjamini-Hochberg over the events tested in the pair) <= 0.05.

Output: ensembl116/events/id50_pairs_fdr/human_mouse_pairs.(tsv|xlsx)
"""
import os

import pandas as pd

import build_supp_tables_116 as bst

E = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50_6ds/'
OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50_pairs_fdr/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
HM = [f'{h}_{m}' for h in ('H1', 'H2', 'H3') for m in ('M1', 'M2')]


def main():
    os.makedirs(OUT, exist_ok=True)
    e = pd.read_csv(E + 'effects/event_effects.tsv', sep='\t',
                    usecols=['cluster', 'comparison', 'cell', 'p.adjust', 'max_abs_dpsi', 'lc_success'])
    e = e[e['comparison'].isin(HM) & e['cell'].isin(CELLS) & e['lc_success'] & e['p.adjust'].notna()]
    rows = []
    for cell in CELLS:
        for pair in HM:
            d = e[(e['cell'] == cell) & (e['comparison'] == pair)]
            sig = d[d['p.adjust'] <= 0.05]
            sd = sig[sig['max_abs_dpsi'] >= 0.1]
            r = {'cell': cell, 'pair': pair.replace('_', '-'), 'tested': len(d)}
            for lab, x in (('sig_FDR_le_0.05', sig), ('sig_FDR_le_0.05_and_abs_max_dPSI_ge_0.1', sd)):
                r[lab] = len(x)
                r[f'pct_{lab}'] = round(100 * len(x) / len(d), 1)
            rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(OUT + 'human_mouse_pairs.tsv', sep='\t', index=False)
    bst.OUT = OUT
    bst.write('human_mouse_pairs', 'Single human-mouse dataset pairs: significant events',
              'tested: events tested by LeafCutter in the pair. Significant: FDR (LeafCutter p.adjust, '
              'Benjamini-Hochberg over the tested events of the pair) <= 0.05, and in the last columns also '
              'max|dPSI| >= 0.1; pct: of the tested events.',
              {'pairs': t})
    pd.set_option('display.width', 250)
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()
