#!/usr/bin/env python3
"""Mouse reads of the event junctions whose mouse partner the liftover changes: old vs proposed mouse coordinates
(raw mouse junction tables of the immune studies). Output: liftover_check/changed_junction_reads.tsv"""
import sys
import pandas as pd
sys.path.insert(0, '/gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116/six')
import check_boundary_mapping_116 as cb

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
OUT = B + 'liftover_check/'
old = pd.read_csv(B + 'orthology/unique_points_116_id50_6ds.txt', sep='\t', index_col=0)
new = pd.read_csv(OUT + 'unique_points_116_id50_liftover.txt', sep='\t', index_col=0)
o, n = dict(zip(old['position_h'], old['position_m'])), dict(zip(new['position_h'], new['position_m']))
changed = {h for h in o if o[h] != n[h]}
ids = pd.read_csv(B + 'events/id50_6ds_anova/AS_clusters_value.txt', sep='\t', usecols=[0]).iloc[:, 0]
rows = []
for j in ids:
    c, s, e, clu = j.split(':')
    hs, he = f'{c}:{s}', f'{c}:{e}'
    if hs in changed or he in changed:
        def mj(m):
            a, b = sorted([int(m[hs].split(':')[1]), int(m[he].split(':')[1])])
            return f"{m[hs].split(':')[0]}:{a}:{b}"
        rows.append(dict(cluster=clu, junction_human=f'{c}:{s}:{e}', mouse_old=mj(o), mouse_new=mj(n)))
t = pd.DataFrame(rows)
ends = set()
for col in ('mouse_old', 'mouse_new'):
    for x in t[col]:
        c, a, b = x.split(':')
        ends |= {(c, int(a)), (c, int(b))}
raw = cb.raw_near('mouse', ends)
t['reads_old'] = t['mouse_old'].map(raw).fillna(0).astype(int)
t['reads_new'] = t['mouse_new'].map(raw).fillna(0).astype(int)
t.to_csv(OUT + 'changed_junction_reads.tsv', sep='\t', index=False)
print('junctions:', len(t), 'events:', t['cluster'].nunique())
print('more reads with new:', (t.reads_new > t.reads_old).sum(), '| fewer:', (t.reads_new < t.reads_old).sum(),
      '| equal:', (t.reads_new == t.reads_old).sum())
print('old ~0 (<10) -> new >=100:', ((t.reads_old < 10) & (t.reads_new >= 100)).sum(),
      '| old >=100 -> new <10:', ((t.reads_old >= 100) & (t.reads_new < 10)).sum())
print(t.sort_values('reads_new', ascending=False).to_string(index=False))
