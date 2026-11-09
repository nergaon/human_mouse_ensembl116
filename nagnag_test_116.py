#!/usr/bin/env python3
"""Are alternative 3'/5' splice site events with small shifts (<= 6 nt, NAGNAG-like) more often species-divergent?
Output: events/id50/revision/small_shift_test.tsv (Fisher's exact test per definition and event type)."""
import numpy as np
import pandas as pd
from scipy import stats

ED = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50/'
s = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
a = pd.read_csv(ED + 'classification/events_annotation.tsv', sep='\t').set_index('cluster')


def shift(r):
    if r.event_type not in ('A3SS', 'A5SS'):
        return np.nan
    (s1, e1), (s2, e2) = [tuple(map(int, j.split('-'))) for j in r.junctions.split(';')]
    return abs(e1 - e2) if s1 == s2 else abs(s1 - s2)


s = s.join(a.apply(shift, axis=1).rename('shift_nt'), on='cluster')
t = s[s['status_high'] != 'not testable']
rows = []
for et in ('A3SS', 'A5SS'):
    x = t[t['event_type'] == et]
    small = x['shift_nt'] <= 6
    for n in ('high', 'high_low'):
        d = x[f'status_{n}'] == 'divergent'
        odds, p = stats.fisher_exact(pd.crosstab(small, d).values)
        rows.append(dict(event_type=et, definition=n, small_shift_event_cells=int(small.sum()),
                         pct_divergent_small=100 * d[small].mean(), larger_shift_event_cells=int((~small).sum()),
                         pct_divergent_larger=100 * d[~small].mean(), odds_ratio=odds, fisher_p=p))
r = pd.DataFrame(rows).round(4)
r.to_csv(ED + 'revision/small_shift_test.tsv', sep='\t', index=False)
print(r.to_string(index=False))
