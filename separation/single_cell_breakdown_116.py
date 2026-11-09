#!/usr/bin/env python3
"""Status, in the other immune cell types, of events that are species-divergent with high confidence in only
one cell type (is 'high confidence in one cell type' evidence of cell type specificity?).
Output: events/id50_sep/revision/single_cell_type_divergent_breakdown.tsv"""
import pandas as pd

ED = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50_sep/'
CROSS = ['H1_M1', 'H1_M2', 'H2_M1', 'H2_M2', 'H3_M1', 'H3_M2']
s = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')


def pairs_summary(x):
    if (x == 'not_testable').all():
        return 'not testable in any human-mouse pair'
    if (x == 'equivalent').any():
        return 'equivalent in >= 1 human-mouse pair'
    if (x[x != 'not_testable'] == 'DS').all():
        return 'differentially spliced in all testable human-mouse pairs'
    if (x == 'DS').any():
        return 'differentially spliced in some human-mouse pairs'
    return 'inconclusive in all testable human-mouse pairs'


rows = []
for name in ('high', 'high_low'):
    col = f'status_{name}'
    div = s[s[col] == 'divergent']
    n = div.groupby('cluster')['cell'].nunique()
    single = div[div['cluster'].isin(n[n == 1].index)][['cluster', 'cell']].rename(columns={'cell': 'hc_cell'})
    other = s.merge(single, on='cluster')
    other = other[other['cell'] != other['hc_cell']]
    e = ev[ev['comparison'].isin(CROSS)].merge(single, on='cluster')
    e = e[e['cell'] != e['hc_cell']]
    pairs = e.groupby(['cluster', 'cell'])['status_m0.1'].agg(pairs_summary)
    base = dict(definition=name, divergent_events=len(n), events_high_confidence_in_one_cell_type=len(single),
                event_x_other_cell_type=len(other))
    for k, v in other['status_high'].value_counts(normalize=True).items():
        rows.append({**base, 'measure': 'high-confidence status in other cell types', 'category': k, 'percent': 100 * v})
    for k, v in other['status_high_low'].value_counts(normalize=True).items():
        rows.append({**base, 'measure': 'high+low-confidence status in other cell types', 'category': k, 'percent': 100 * v})
    for k, v in pairs.value_counts(normalize=True).items():
        rows.append({**base, 'measure': 'human-mouse pairs in other cell types', 'category': k, 'percent': 100 * v})
r = pd.DataFrame(rows).round(1)
r.to_csv(ED + 'revision/single_cell_type_divergent_breakdown.tsv', sep='\t', index=False)
print(r[r.definition == 'high'].to_string(index=False))
