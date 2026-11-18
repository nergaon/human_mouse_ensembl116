#!/usr/bin/env python3
"""
Exon labels of the boundary pairs switched by the liftOver rule (liftover_boundaries.py).

A switched pair keeps its human boundary and gets a new mouse boundary, but the exon labels and scores of the row were
those of the old mouse partner. For each switched row (orthology/unique_points_116_<variant>.txt vs the backup in
pre_liftover_20261118/orthology/):
  1. an exon pair of the variant (exon_pairs_116.tsv, same filters as orthology_points_116.py) with the same human
     boundary and the new mouse boundary: labels and scores from the most canonical such pair (as select_best);
  2. otherwise the mouse exon with that boundary (exons_m_116.tsv, the mouse ortholog gene; canonical rank first):
     exon_m and rank_m only (the scores stay those of the human exon's best pair).
Writes the corrected files in place (atomic rename) and liftover_check/<variant>/label_fixes.tsv.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, '/gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116')
from orthology_points_116 import POINT_COLUMNS, VARIANTS, to_points

E = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
O = E + 'orthology/'
LABEL = ['exon_h', 'rank_h', 'exon_m', 'rank_m', '%id_human', '%id_mouse', 'Average %id', 'Ratio_based_score']


def main():
    pairs = pd.read_csv(O + 'exon_pairs_116.tsv', sep='\t')
    pairs['min_id'] = pairs[['%id_human', '%id_mouse']].min(axis=1)
    exons_m = pd.read_csv(O + 'exons_m_116.tsv', sep='\t')
    gene_m = pairs.drop_duplicates('ensembl').set_index('ensembl')['ensembl_m']
    chrom_m = pairs.assign(c=pairs['start_m'].str.split(':').str[0]).drop_duplicates('ensembl_m') \
        .set_index('ensembl_m')['c']
    for v, (min_id, rbh) in VARIANTS.items():
        cur = pd.read_csv(O + f'unique_points_116_{v}.txt', sep='\t', index_col=0)
        old = pd.read_csv(E + f'pre_liftover_20261118/orthology/unique_points_116_{v}.txt', sep='\t', index_col=0)
        sw = cur.index[cur['position_m'] != old['position_m']]
        p = pairs[pairs['min_id'] >= min_id]
        if rbh:
            p = p[p['reciprocal_best']]
        pts = to_points(p)
        pts['star_count'] = pts['rank_h'].astype(str).add(pts['rank_m'].astype(str)).str.count(r'\*')
        pts = pts.sort_values(['star_count', 'Ratio_based_score'], ascending=[True, False]) \
            .drop_duplicates(['position_h', 'position_m']).set_index(['position_h', 'position_m'])
        rows = []
        for i in sw:
            r = cur.loc[i]
            key = (r['position_h'], r['position_m'])
            if key in pts.index:
                x = pts.loc[key]
                for c in LABEL + ['star_count']:
                    cur.at[i, c] = x[c]
                how = 'exon pair'
            else:
                gm = gene_m.get(r['ensembl'])
                c, pos = r['position_m'].split(':')
                em = exons_m[(exons_m['ensembl'] == gm) & ((exons_m['start'] == int(pos)) | (exons_m['end'] == int(pos)))]
                if gm is None or em.empty or chrom_m.get(gm) != c:
                    rows.append(dict(row=i, position_h=r['position_h'], position_m=r['position_m'], how='not found'))
                    continue
                em = em.assign(star=em['rank'].astype(str).str.count(r'\*')).sort_values('star').iloc[0]
                cur.at[i, 'exon_m'], cur.at[i, 'rank_m'] = em['exon'], em['rank']
                cur.at[i, 'star_count'] = str(r['rank_h']).count('*') + str(em['rank']).count('*')
                how = 'mouse exon'
            rows.append(dict(row=i, position_h=r['position_h'], position_m=r['position_m'], how=how,
                             old_exon_m=old.at[i, 'exon_m'], old_rank_m=old.at[i, 'rank_m'],
                             exon_m=cur.at[i, 'exon_m'], rank_m=cur.at[i, 'rank_m']))
        t = pd.DataFrame(rows)
        d = E + 'liftover_check/' + ('' if v == 'id50' else v + '/')
        t.to_csv(d + 'label_fixes.tsv', sep='\t', index=False)
        tmp = O + f'.unique_points_116_{v}.tmp'
        cur[POINT_COLUMNS].to_csv(tmp, sep='\t')
        os.replace(tmp, O + f'unique_points_116_{v}.txt')
        print(v, len(sw), 'switched rows:', t['how'].value_counts().to_dict(),
              '| rank_m no longer *:', int(((t.get('old_rank_m') == '*') & (t.get('rank_m') != '*')).sum()))


if __name__ == '__main__':
    main()
