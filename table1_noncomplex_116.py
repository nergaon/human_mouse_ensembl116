#!/usr/bin/env python3
"""
Table 1 of the revised manuscript restricted to the non-complex events (SE, A3SS, A5SS), as a sensitivity
check; not in the manuscript. Same rows and definitions as table1_rows() of make_revised_docx_116.py
(classification counts from event_status_by_cell.tsv, per-pair rates as in facts_116.py).
Output: ensembl116/sensitivity/Table1_noncomplex.docx / .tsv (and Table1_all.tsv, the same code on all events,
which must reproduce the manuscript Table 1)
"""
import os

import docx
import pandas as pd
from docx.shared import Pt

from docx_style_116 import finalize

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50/'
OUT = B + 'sensitivity/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
DS_LABEL = 'Differentially spliced events [p≤0.05 & max(|ΔPSI|)≥0.1], %'
EQ_LABEL = 'Equivalent events [90% CI of ΔPSI within ±0.1], %'


def kind(c):  # as in facts_116.py (not imported: that module writes its outputs on import)
    if c.startswith('matched_'):
        return 'matched_' + kind(c[8:])
    if c == 'FH_FM':
        return 'fibro'
    a, b = c.split('_')
    return 'cross' if a[0] != b[0] else ('HH' if a[0] == 'H' else 'MM')


def pair_rates(t, k, cell):
    g = t[(t['kind'] == k) & (t['cell'] == cell)]
    if g.empty:
        return None
    per = g.groupby('comparison').agg(n=('cluster', 'size'),
                                      DS=('status_m0.1', lambda s: 100 * (s == 'DS').mean()),
                                      EQ=('status_m0.1', lambda s: 100 * (s == 'equivalent').mean()),
                                      INC=('status_m0.1', lambda s: 100 * (s == 'inconclusive').mean()))
    return {c: (f'{per[c].median():.0f} [{per[c].min():.0f}-{per[c].max():.0f}]' if len(per) > 1
                else f'{per[c].iloc[0]:.0f}') for c in per}


def rows(clusters):
    st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
    st = st[st['cluster'].isin(clusters)]
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    ev = ev[ev['cluster'].isin(clusters)]
    ev['kind'] = ev['comparison'].map(kind)
    t = ev[ev['status_m0.1'] != 'not_testable']
    test = st[st['status_high'] != 'not testable']
    tc = test.groupby('cell').size()

    def hcrow(col, s):
        n = test[test[col] == s].groupby('cell').size()
        return [f'{n.get(c, 0)} ({100 * n.get(c, 0) / tc[c]:.1f}%)' for c in CELLS]
    fs = t.loc[t['kind'] == 'fibro', 'status_m0.1'].value_counts()
    fn = int(fs.sum())

    def fcell(k):
        return f'{fs.get(k, 0)} ({100 * fs.get(k, 0) / fn:.1f}%)'

    def prow(label, k, col):
        return ('row', [label] + [(pair_rates(t, k, c) or {}).get(col, '') for c in CELLS] + [''])
    R = [('header', ['Comparison category', 'CD4T', 'CD8T', 'B', 'NK', 'Mono', 'Fibroblasts']),
         ('section', ['Human-Mouse alternative splicing events (immune cells: 6 dataset pairs, classified per event '
                      'from two pairs that share no dataset; fibroblasts: 1 pair)']),
         ('row', ['Testable events'] + [str(tc[c]) for c in CELLS] + [str(fn)]),
         ('row', ['Differentially spliced events, high confidence'] + hcrow('status_high', 'divergent') + ['–']),
         ('row', ['Differentially spliced events, high + low confidence'] + hcrow('status_high_low', 'divergent')
          + [fcell('DS')]),
         ('row', ['Equivalent events, high confidence'] + hcrow('status_high', 'conserved') + ['–']),
         ('row', ['Equivalent events, high + low confidence'] + hcrow('status_high_low', 'conserved')
          + [fcell('equivalent')]),
         ('row', ['Inconclusive events'] + hcrow('status_high_low', 'inconclusive') + [fcell('inconclusive')]),
         ('section', ['Human-Mouse, matched populations (M1 naive T cells, blood monocytes; 3 pairs)']),
         prow('Differentially spliced events, %', 'matched_cross', 'DS')]
    for k, title in (('HH', 'Human-Human alternative splicing events (3 dataset pairs)'),
                     ('MM', 'Mouse-Mouse alternative splicing events (1 dataset pair)')):
        R += [('section', [title]), prow('Comparable events', k, 'n'), prow(DS_LABEL, k, 'DS'),
              prow(EQ_LABEL, k, 'EQ'), prow('Inconclusive events, %', k, 'INC')]
    return R


def to_tsv(R, path):
    pd.DataFrame([v + [''] * (7 - len(v)) for _, v in R]).to_csv(path, sep='\t', index=False, header=False)


def to_docx(R, caption, path):
    d = docx.Document()
    d.add_paragraph().add_run(caption).bold = True
    tab = d.add_table(rows=0, cols=7)
    tab.style = 'Table Grid'
    for k, vals in R:
        cells = tab.add_row().cells
        if k == 'section':
            cells[0].merge(cells[6]).text = vals[0]
        else:
            for c, v in zip(cells, vals):
                c.text = str(v)
        for c in cells:
            for p in c.paragraphs:
                p.paragraph_format.space_after = Pt(0)
                for r in p.runs:
                    r.bold = k in ('header', 'section')
    finalize(d)
    d.save(path)


def main():
    os.makedirs(OUT, exist_ok=True)
    ann = pd.read_csv(ED + 'classification/events_annotation.tsv', sep='\t')
    simple = ann.loc[ann['event_type'] != 'complex', 'cluster']
    to_tsv(rows(ann['cluster']), OUT + 'Table1_all.tsv')
    R = rows(simple)
    to_tsv(R, OUT + 'Table1_noncomplex.tsv')
    to_docx(R, f'Table 1, non-complex events only ({len(simple):,} of {len(ann):,} orthologous alternative splicing '
               'events: SE, A3SS and A5SS; complex events excluded). Same definitions as Table 1.',
            OUT + 'Table1_noncomplex.docx')


if __name__ == '__main__':
    main()
