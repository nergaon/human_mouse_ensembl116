#!/usr/bin/env python3
"""
Orthologous exon boundaries ("points") from exon_pairs_116.tsv, for several
orthology stringencies (reviewer request: sensitivity to the 50% identity cutoff).

Replaces sum_data_HN12.py. For each variant, pairs are filtered, genes with fewer
than 2 orthologous exon pairs are dropped, every pair contributes its 5' and 3'
boundaries, and select_best keeps a human-mouse boundary pair only if each boundary is the
other's best partner (mutual best: highest Ratio_based_score, then fewer non-canonical exons).
(Before 2026-11-02 the best partner was chosen per human boundary and then per mouse boundary
by score only; a mouse boundary whose best human partner had been given another mouse
boundary was then mapped to a worse, often non-orthologous, human boundary - e.g. the
canonical FERMT3 E12 acceptor - creating false species-specific splice sites.)

Variants:
  id50, id70, id80   min(%id_human, %id_mouse) >= 50 / 70 / 80
  id50_rbh           id50 and the pair is reciprocal best
Outputs (ensembl116/orthology/):
  unique_points_116_<variant>.txt   same columns as ensembl115/unique_points_HN6.txt
  orthology_summary_116.tsv         gene / exon / boundary counts per variant
"""
import pandas as pd

ORTH_DIR = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/orthology/'
VARIANTS = {'id50': (50, False), 'id70': (70, False), 'id80': (80, False), 'id50_rbh': (50, True)}
POINT_COLUMNS = ['position_m', 'position_h', 'ensembl', 'symbol', 'exon_h', 'rank_h', 'exon_m', 'rank_m',
                 '%id_human', '%id_mouse', 'Average %id', 'Ratio_based_score', 'star_count']


def select_best(df):
    df = df.copy()
    df['star_count'] = df['rank_h'].astype(str).add(df['rank_m'].astype(str)).str.count(r'\*')
    df = df.sort_values(['Ratio_based_score', 'star_count', 'Average %id', 'position_h', 'position_m'],
                        ascending=[False, True, False, True, True])
    best_h = df.drop_duplicates('position_h')                     # best mouse partner of each human boundary
    best_m = df.drop_duplicates('position_m')[['position_h', 'position_m']]   # best human partner of each mouse one
    best = best_h.merge(best_m, on=['position_h', 'position_m'])
    # exon labels (exon_h/m, rank_h/m) from the most canonical exon pair with these same two boundaries (pairing
    # unchanged): the best-scoring pair is often a short non-canonical exon sharing the boundary, e.g. HSP90B1 E3
    # donor labeled '*' from a 28-nt CDS fragment (score 96) instead of E3 itself (92) (fixed 2026-11-03)
    label = ['exon_h', 'rank_h', 'exon_m', 'rank_m', 'star_count']
    canon = df.sort_values(['star_count', 'Ratio_based_score'], ascending=[True, False]) \
        .drop_duplicates(['position_h', 'position_m'])[['position_h', 'position_m'] + label]
    best = best.drop(columns=label).merge(canon, on=['position_h', 'position_m'])
    return best.sort_values('position_m').reset_index(drop=True)


def to_points(pairs):
    common = ['ensembl', 'symbol', 'exon_h', 'rank_h', 'exon_m', 'rank_m',
              '%id_human', '%id_mouse', 'Average %id', 'Ratio_based_score']
    starts = pairs[common + ['start_h', 'start_m']].rename(columns={'start_h': 'position_h', 'start_m': 'position_m'})
    ends = pairs[common + ['end_h', 'end_m']].rename(columns={'end_h': 'position_h', 'end_m': 'position_m'})
    return pd.concat([starts, ends], ignore_index=True)


def main():
    pairs = pd.read_csv(ORTH_DIR + 'exon_pairs_116.tsv', sep='\t')
    exons_h = pd.read_csv(ORTH_DIR + 'exons_h_116.tsv', sep='\t')
    exons_m = pd.read_csv(ORTH_DIR + 'exons_m_116.tsv', sep='\t')
    status = pd.read_csv(ORTH_DIR + 'genes_status_116.tsv', sep='\t')
    multi_exon = status[(status['exons_h'] >= 2) & (status['exons_m'] >= 2)]
    pairs['min_id'] = pairs[['%id_human', '%id_mouse']].min(axis=1)

    summary = []
    for name, (min_id, rbh) in VARIANTS.items():
        p = pairs[pairs['min_id'] >= min_id]
        if rbh:
            p = p[p['reciprocal_best']]
        n_pairs = p.groupby('ensembl').size()
        p = p[p['ensembl'].isin(n_pairs[n_pairs > 1].index)]
        points = select_best(to_points(p))[POINT_COLUMNS]
        points.to_csv(ORTH_DIR + f'unique_points_116_{name}.txt', sep='\t')
        genes = set(p['ensembl'])
        summary.append({
            'variant': name,
            'multi_exon_ortholog_genes': len(multi_exon),
            'genes_with_ortholog_exons': len(genes),
            'exons_h': int(exons_h['ensembl'].isin(genes).sum()),
            'exons_h_with_ortholog': p['exon_h'].nunique(),
            'exons_m': int(exons_m['ensembl'].isin(multi_exon['ensembl_m'][multi_exon['ensembl_h'].isin(genes)]).sum()),
            'exons_m_with_ortholog': p['exon_m'].nunique(),
            'exon_pairs': len(p),
            'unique_boundaries': len(points),
        })
        print(summary[-1], flush=True)
    pd.DataFrame(summary).to_csv(ORTH_DIR + 'orthology_summary_116.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
