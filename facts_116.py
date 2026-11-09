#!/usr/bin/env python3
"""
All numbers quoted in the revised manuscript (Ensembl 116, id50 unless stated), written to
ensembl116/revision_files/manuscript_numbers.md, plus the new Table 1 (table1.tsv).
"""
import numpy as np
import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
out = []


def p(*a):
    out.append(' '.join(str(x) for x in a))


def kind(c):
    if c.startswith('matched_'):
        return 'matched_' + kind(c[8:])
    if c == 'FH_FM':
        return 'fibro'
    a, b = c.split('_')
    return 'cross' if a[0] != b[0] else ('HH' if a[0] == 'H' else 'MM')


orth = pd.read_csv(B + 'orthology/orthology_summary_116.tsv', sep='\t')
st = pd.read_csv(B + 'orthology/genes_status_116.tsv', sep='\t')
p('## Orthology'); p(orth.to_string(index=False))
p('ortholog gene pairs:', len(st), st['status'].value_counts().to_dict())
ann = pd.read_csv(ED + 'classification/events_annotation.tsv', sep='\t')
p('## Events'); p('events', len(ann), 'genes', ann['ensembl'].nunique(), 'junctions', ann['n_junctions'].sum(),
                  ann['event_type'].value_counts().to_dict())
p('splice-site correspondence', ann['splice_site_correspondence'].value_counts(dropna=False).to_dict())
pres = pd.read_csv(ED + 'junction_presence.tsv', sep='\t', index_col=0)
p('orthologous junctions observed (immune)', len(pres), 'kept', int(pres['kept'].sum()))
p('genes with 1 event', int((ann.groupby('ensembl').size() == 1).sum()), 'with >1', int((ann.groupby('ensembl').size() > 1).sum()),
  'max events per gene', int(ann.groupby('ensembl').size().max()))

ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
ev['kind'] = ev['comparison'].map(kind)
status = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
testable_cells = status[status['status_high'] != 'not testable']
p('## Testable'); p('events testable in two independent human-mouse pairs in >=1 cell type:', testable_cells['cluster'].nunique(),
                    'genes', testable_cells['ensembl'].nunique())
p('per cell:', testable_cells.groupby('cell').size().to_dict())

# Table 1: median [min-max] over pairs
rows = []
t = ev[ev['status_m0.1'] != 'not_testable']
for k in ('cross', 'HH', 'MM', 'fibro', 'matched_cross'):
    for cell in CELLS + ['Fibro']:
        g = t[(t['kind'] == k) & (t['cell'] == cell)]
        if g.empty:
            continue
        per = g.groupby('comparison').agg(n=('cluster', 'size'),
                                          DS=('status_m0.1', lambda s: 100 * (s == 'DS').mean()),
                                          EQ=('status_m0.1', lambda s: 100 * (s == 'equivalent').mean()),
                                          INC=('status_m0.1', lambda s: 100 * (s == 'inconclusive').mean()),
                                          EQ05=('status_m0.05', lambda s: 100 * (s == 'equivalent').mean()))
        r = dict(kind=k, cell=cell, pairs=len(per))
        for c in per:
            r[c] = (f'{per[c].median():.0f} [{per[c].min():.0f}-{per[c].max():.0f}]' if len(per) > 1
                    else f'{per[c].iloc[0]:.0f}')
        rows.append(r)
tab = pd.DataFrame(rows)
p('## Table 1 (median [range] over dataset pairs; % of comparable events)'); p(tab.to_string(index=False))
hc = pd.read_csv(ED + 'classification/summary_by_cell.tsv', sep='\t', index_col=0)
tc = testable_cells.groupby('cell').size()
for col in ['high_divergent', 'high_conserved', 'high_low_divergent', 'high_low_conserved']:
    p(col, {c: f'{hc.loc[c, col]} ({100 * hc.loc[c, col] / tc[c]:.1f}%)' for c in CELLS})
for n in ('high', 'high_low'):
    for s in ('divergent', 'conserved'):
        sub = status[status[f'status_{n}'] == s]
        cnt = sub.groupby('cluster')['cell'].nunique()
        p(f'{n} {s}: events {len(cnt)}, genes {sub["ensembl"].nunique()}, in 1 cell {int((cnt == 1).sum())}, '
          f'in >1 {int((cnt > 1).sum())}, in all 5 {int((cnt == 5).sum())}:',
          ', '.join(sorted(sub[sub['cluster'].isin(cnt[cnt == 5].index)]['symbol'].unique())))
        p('   unique to one cell type by cell:', sub[sub['cluster'].isin(cnt[cnt == 1].index)].groupby('cell').size().to_dict())
tab.to_csv(B + 'revision_files/table1.tsv', sep='\t', index=False)

# comparison type summaries
s = pd.read_csv(ED + 'effects/comparison_type_summary.tsv', sep='\t')
p('## Comparison types'); p(s.round(2).to_string(index=False))
for v in ['id70', 'id80', 'id50_rbh']:
    c = pd.read_csv(B + f'events/{v}/classification/summary_by_cell.tsv', sep='\t', index_col=0)
    a = pd.read_csv(B + f'events/{v}/classification/events_annotation.tsv', sep='\t')
    st_v = pd.read_csv(B + f'events/{v}/classification/event_status_by_cell.tsv', sep='\t')
    ts = pd.read_csv(B + f'events/{v}/effects/comparison_type_summary.tsv', sep='\t')
    cr = ts[ts['kind'] == 'cross']
    p(f'{v}: events {len(a)}, high div events {st_v.loc[st_v.status_high == "divergent", "cluster"].nunique()}, '
      f'high cons events {st_v.loc[st_v.status_high == "conserved", "cluster"].nunique()}, '
      f'cross DS {cr.pct_DS_min.min():.1f}-{cr.pct_DS_max.max():.1f}, cross r {cr.pearson_min.min():.2f}-{cr.pearson_max.max():.2f}')
for f in ['dpsi_cutoff_sensitivity', 'common_events_tests', 'splice_site_subset_summary', 'coverage_bias_model',
          'fibroblast_effects', 'enrichment_summary', 'minor_usage_summary']:
    p(f'## {f}'); p(pd.read_csv(ED + f'revision/{f}.tsv', sep='\t').round(4).to_string(index=False))
ce = pd.read_csv(ED + 'revision/common_events_summary.tsv', sep='\t')
p('## common events by kind'); p(ce.groupby('kind')[['common_events', 'pct_DS', 'pct_equivalent_m005', 'median_max_abs_dpsi', 'pearson']]
                                 .agg(['min', 'max']).round(2).to_string())
# fibroblasts (Figure 6 A-D numbers)
hm = ev[ev['comparison'] == 'H2_M1']
fib = ev[ev['comparison'] == 'FH_FM']
ft = fib[fib['status_m0.1'] != 'not_testable']
p('## Fibroblasts'); p('comparable', len(ft), (ft['status_m0.1'].value_counts(normalize=True) * 100).round(1).to_dict(),
                       ft['status_m0.1'].value_counts().to_dict())
fmap = fib.set_index('cluster')['status_m0.1']
cls = status[status['status_high_low'] != 'not testable']
for k in ('divergent', 'conserved'):
    for cell in CELLS:
        b = cls[(cls.cell == cell) & (cls.status_high_low == k)]['cluster'].map(fmap)
        b = b[b.notna() & (b != 'not_testable')]
        p(f'immune {k} {cell} -> fibro', len(b), (b.value_counts(normalize=True) * 100).round(0).to_dict())
for base in ('DS', 'equivalent'):
    fc = set(fib[fib['status_m0.1'] == base]['cluster'])
    for cell in CELLS:
        im = cls[(cls.cell == cell) & cls.cluster.isin(fc)]['status_high_low']
        p(f'fibro {base} -> {cell}', len(im), (im.value_counts(normalize=True) * 100).round(0).to_dict())
for n in ('high', 'high_low'):
    for sname in ('divergent', 'conserved'):
        cl = status.loc[status[f'status_{n}'] == sname, 'cluster'].unique()
        f = pd.Series(cl).map(fmap)
        f = f[f.notna() & (f != 'not_testable')]
        p(f'{n} {sname} -> fibro', len(f), (f.value_counts(normalize=True) * 100).round(1).to_dict())
open(B + 'revision_files/manuscript_numbers.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
