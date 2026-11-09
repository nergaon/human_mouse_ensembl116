#!/usr/bin/env python3
"""Exon rank consistency of the orthologous exon boundaries used by the AS events (rank is not used to define
orthology; this is a validation). Rank = exon number in the canonical transcript of each species.
Output: orthology/rank_consistency_116_<variant>.tsv"""
import re
import sys

import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
v = sys.argv[1] if len(sys.argv) > 1 else 'id50'
pts = pd.read_csv(B + f'orthology/unique_points_116_{v}.txt', sep='\t', index_col=0)
ann = pd.read_csv(B + f'events/{v}/classification/events_annotation.tsv', sep='\t')
rows = [(r.ensembl, f'{r.chrom}:{p}') for r in ann.itertuples() for j in r.junctions.split(';') for p in j.split('-')]
b = pd.DataFrame(rows, columns=['ensembl', 'position_h']).drop_duplicates('position_h')
b = b.merge(pts[['position_h', 'rank_h', 'rank_m']], on='position_h', how='left')
num = lambda x: int(m.group(1)) if (m := re.match(r'E(\d+)', str(x))) else None
b['rh'], b['rm'] = b['rank_h'].map(num), b['rank_m'].map(num)
both = b.dropna(subset=['rh', 'rm']).copy()
both['offset'] = both['rh'] - both['rm']
mode = both.groupby('ensembl')['offset'].agg(lambda x: x.mode().iloc[0])
res = pd.DataFrame([dict(
    orthology=v, event_boundaries=len(b), in_canonical_transcript_of_both=len(both),
    pct_in_canonical_of_both=100 * len(both) / len(b), pct_same_rank=100 * (both['offset'] == 0).mean(),
    pct_constant_offset_within_gene=100 * (both['offset'] == both['ensembl'].map(mode)).mean(),
    pct_rank_difference_le_1=100 * (both['offset'].abs() <= 1).mean())]).round(1)
res.to_csv(B + f'orthology/rank_consistency_116_{v}.tsv', sep='\t', index=False)
print(res.to_string(index=False))
