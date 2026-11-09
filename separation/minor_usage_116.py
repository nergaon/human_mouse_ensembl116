#!/usr/bin/env python3
"""
Minor-isoform usage of each event x cell type and the classification stratified by it.

usage = PSI of the less used isoform, averaged over the human and mouse groups of the six
human-mouse comparisons:
  SE       exclusion = skip / (skip + mean(two inclusion junctions)); usage = min(incl, excl)
  A3/A5SS  smaller of the two junction fractions
  complex  1 - largest junction fraction
Usage: minor_usage_116.py <variant>  -> events/<variant>/revision/minor_usage*.tsv
"""
import sys

import numpy as np
import pandas as pd

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/'
CROSS = ['H1_M1', 'H1_M2', 'H2_M1', 'H2_M2', 'H3_M1', 'H3_M2']


def usage(group, etype, junctions):
    p = group.set_index('intron')
    vals = []
    for col in ('psi1', 'psi2'):
        x = p[col]
        if x.isna().all():
            continue
        if etype == 'SE':
            ends = {i: tuple(map(int, i.split(':')[1:3])) for i in x.index}
            skip = max(ends, key=lambda i: ends[i][1] - ends[i][0])
            incl = x.drop(skip).mean()
            excl = x[skip] / (x[skip] + incl) if (x[skip] + incl) > 0 else np.nan
            vals.append(min(excl, 1 - excl))
        elif etype in ('A3SS', 'A5SS'):
            vals.append(x.min())
        else:
            vals.append(1 - x.max())
    return np.nanmean(vals) if vals else np.nan


def main():
    v = sys.argv[1]
    s = pd.read_csv(BASE + f'{v}/classification/event_status_by_cell.tsv', sep='\t')
    j = pd.read_csv(BASE + f'{v}/effects/junction_effects.tsv', sep='\t')
    j = j[j['comparison'].isin(CROSS) & j['testable']]
    etype = s.drop_duplicates('cluster').set_index('cluster')['event_type']
    u = j.groupby(['cell', 'cluster', 'comparison']).apply(lambda g: usage(g, etype[g.name[1]], None))
    u = u.groupby(['cell', 'cluster']).mean().rename('minor_usage')
    s = s.join(u, on=['cell', 'cluster'])
    s.to_csv(BASE + f'{v}/revision/minor_usage_by_event.tsv', sep='\t', index=False)
    t = s[s['status_high'] != 'not testable'].copy()
    t['usage_bin'] = pd.cut(t['minor_usage'], [-0.001, 0.01, 0.05, 0.1, 0.2, 0.5],
                            labels=['<1%', '1-5%', '5-10%', '10-20%', '20-50%'])
    out = []
    for col in ('status_high', 'status_high_low'):
        n = pd.crosstab(t['usage_bin'], t[col])
        pct = pd.crosstab(t['usage_bin'], t[col], normalize='index').mul(100).round(1).add_suffix(' %')
        out.append(pd.concat([n, pct], axis=1).assign(definition=col.replace('status_', '')))
    res = pd.concat(out).reset_index()
    res.to_csv(BASE + f'{v}/revision/minor_usage_summary.tsv', sep='\t', index=False)
    print(res.to_string(index=False))


if __name__ == '__main__':
    main()
