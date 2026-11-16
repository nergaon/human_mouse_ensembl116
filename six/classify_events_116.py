#!/usr/bin/env python3
"""
Event annotation and human-mouse classification with confidence levels (revision).

Per AS event (LeafCutter cluster): gene, event type, human and orthologous mouse coordinates, and whether the
splice-site correspondence is well supported (A3/A5: identical distance between the alternative sites in
both species; skipped exon: identical exon length).
Event types (human Ensembl 116 annotation, all transcripts of the gene):
  SE        three junctions: inclusion, inclusion, skipping
  A3SS/A5SS two junctions sharing one end, alternative sites on the same exon (an annotated exon contains both)
  AFE/ALE   two junctions sharing one end, alternative sites on different exons, one of which is a first (AFE) or
            last (ALE) exon of an annotated transcript (alternative first / last exon)
  complex   all other events (including two-junction events with alternative internal exons)

Per cell type, from the comparison statuses of effect_size_analysis.py
(DS = p.adjust <= 0.05 and max|dPSI| >= 0.1; equivalent = every junction's 90% CI of dPSI within +-0.1):
  testable   at least two human-mouse pairs that share no dataset (e.g. H1-M1 and H2-M2) are testable
  divergent  DS in two human-mouse pairs that share no dataset, and equivalent in no human-mouse pair
  conserved  equivalent in two human-mouse pairs that share no dataset, and DS in no human-mouse pair
  confidence (divergent and conserved alike)
    high     no within-species pair DS, and at least one human-human pair and the mouse-mouse pair equivalent
    low      all other divergent / conserved events (the within-species comparisons are not used)
  all other testable events are inconclusive
  status_high      divergent / conserved = high confidence
  status_high_low  divergent / conserved = high or low confidence (nested: includes the high-confidence events)
The manuscript's original M1 groups are used (the 'matched' set is a separate sensitivity analysis).

Usage: classify_events_116.py <variant>     (e.g. id50)
Output (ensembl116/events/<variant>/classification/):
  events_annotation.tsv, event_status_by_cell.tsv, summary_by_cell.tsv
"""
import gzip
import itertools
import math
import os
import re
import sys

import numpy as np
import pandas as pd

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
GTF_H = BASE + 'ref/Homo_sapiens.GRCh38.116.gtf.gz'
CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']
STATUS = 'status_m0.1'  # DS, or equivalent within +-0.1
LEVELS = ('high', 'high_low')


def gene_annotation():
    """Strand per gene and, per gene, its exons (start, end, first exon, last exon) from all transcripts."""
    strands, exons, n_exons = {}, {}, {}
    rows = []
    with gzip.open(GTF_H, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.split('\t', 9)
            if p[2] == 'gene':
                strands[re.search(r'gene_id "([^"]+)"', p[8]).group(1)] = p[6]
            elif p[2] == 'exon':
                g = re.search(r'gene_id "([^"]+)"', p[8]).group(1)
                t = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
                n = int(re.search(r'exon_number "(\d+)"', p[8]).group(1))
                rows.append((g, t, int(p[3]), int(p[4]), n))
                n_exons[t] = max(n_exons.get(t, 0), n)
    for g, t, st, en, n in rows:
        exons.setdefault(g, set()).add((st, en, n == 1, n == n_exons[t]))
    return strands, exons


def event_type(junctions, strand, exons=()):
    """junctions: list of (start, end) with start < end, human coordinates (start = last base of the left exon,
    end = first base of the right exon); exons: (start, end, is_first, is_last) of the gene."""
    if len(junctions) == 2:
        (s1, e1), (s2, e2) = junctions
        if s1 == s2 or e1 == e2:
            left = e1 == e2                          # alternative sites at the ends of the left exons
            x, y = sorted((s1, s2) if left else (e1, e2))
            if not exons or any(st <= x and en >= y for st, en, _, _ in exons):
                if left:
                    return 'A5SS' if strand == '+' else 'A3SS'
                return 'A3SS' if strand == '+' else 'A5SS'
            alt = [ex for ex in exons if (ex[1] if left else ex[0]) in (x, y)]
            upstream = left == (strand == '+')       # the alternative exons are 5' of the shared exon
            if upstream and any(ex[2] for ex in alt):
                return 'AFE'
            if not upstream and any(ex[3] for ex in alt):
                return 'ALE'
            return 'complex'
    if len(junctions) == 3:
        js = sorted(junctions)
        long_ = max(js, key=lambda j: j[1] - j[0])
        rest = [j for j in js if j != long_]
        if len(rest) == 2:
            (a, b), (c, d) = sorted(rest)
            if a == long_[0] and d == long_[1] and b < c:
                return 'SE'
    return 'complex'


def correspondence(etype, junctions, h2m):
    """Well-supported splice-site correspondence between human and the orthologous mouse event."""
    def m(pos):
        return h2m.get(pos)
    if etype in ('A3SS', 'A5SS'):
        (s1, e1), (s2, e2) = junctions
        a, b = ((e1, e2) if s1 == s2 else (s1, s2))
        ma, mb = m(a), m(b)
        if ma is None or mb is None:
            return None
        return abs(a - b) == abs(ma - mb)
    if etype == 'SE':
        js = sorted(junctions)
        long_ = max(js, key=lambda j: j[1] - j[0])
        (_, b), (c, _) = sorted(j for j in js if j != long_)
        mb, mc = m(b), m(c)
        if mb is None or mc is None:
            return None
        return abs(c - b) == abs(mc - mb)
    return None


def annotate(variant):
    counts = pd.read_csv(BASE + f'events/{variant}/AS_clusters_value.txt', sep='\t', index_col=0, usecols=[0])
    points = pd.read_csv(BASE + f'orthology/unique_points_116_{variant}.txt', sep='\t', index_col=0)
    by_pos = points.set_index('position_h')
    h2m = {int(p.split(':')[1]): int(q.split(':')[1]) for p, q in zip(points['position_h'], points['position_m'])}
    strands, gene_exons = gene_annotation()
    rows = []
    jun = pd.Series(counts.index).str.split(':', expand=True)
    jun.columns = ['chrom', 'start', 'end', 'cluster']
    for cluster, g in jun.groupby('cluster'):
        chrom = g['chrom'].iloc[0]
        js = [(int(s), int(e)) for s, e in zip(g['start'], g['end'])]
        pos = [f'{chrom}:{p}' for j in js for p in j]
        hit = by_pos.reindex(pos).dropna(subset=['ensembl'])
        ensembl = hit['ensembl'].mode().iloc[0] if len(hit) else ''
        symbol = hit.loc[hit['ensembl'] == ensembl, 'symbol'].iloc[0] if len(hit) else ''
        strand = strands.get(ensembl, '+')
        et = event_type(js, strand, gene_exons.get(ensembl, set()))
        rows.append(dict(cluster=cluster, chrom=chrom, start=min(j[0] for j in js), end=max(j[1] for j in js),
                         strand=strand, ensembl=ensembl, symbol=symbol, n_junctions=len(js), event_type=et,
                         junctions=';'.join(f'{s}-{e}' for s, e in js),
                         mouse_positions=';'.join(f'{h2m.get(s, "NA")}-{h2m.get(e, "NA")}' for s, e in js),
                         splice_site_correspondence=correspondence(et, js, h2m)))
    return pd.DataFrame(rows)


def classify(events):
    e = events[~events['comparison'].str.startswith('matched_') & (events['comparison'] != 'FH_FM')]
    cross = sorted(c for c in e['comparison'].unique() if c[0] != c.split('_')[1][0])
    hh = sorted(c for c in e['comparison'].unique() if c not in cross and c[0] == 'H')
    mm = sorted(c for c in e['comparison'].unique() if c not in cross and c[0] == 'M')
    within = hh + mm
    st = e.pivot_table(index=['cell', 'cluster'], columns='comparison', values=STATUS, aggfunc='first') \
        .reindex(columns=cross + within).fillna('not_testable')
    dp = e.pivot_table(index=['cell', 'cluster'], columns='comparison', values='max_abs_dpsi', aggfunc='first')
    ds, eq, tst = st == 'DS', st == 'equivalent', st != 'not_testable'
    independent = [(a, b) for a, b in itertools.combinations(cross, 2)
                   if a.split('_')[0] != b.split('_')[0] and a.split('_')[1] != b.split('_')[1]]

    def replicated(mask):
        return np.logical_or.reduce([mask[a] & mask[b] for a, b in independent])

    out = pd.DataFrame(index=st.index)
    out['cross_testable'] = tst[cross].sum(1)
    out['cross_DS'] = ds[cross].sum(1)
    out['cross_equivalent'] = eq[cross].sum(1)
    out['within_testable'] = tst[within].sum(1)
    out['within_DS'] = ds[within].sum(1)
    out['within_equivalent'] = eq[within].sum(1)
    out['median_cross_abs_dpsi'] = dp[cross].where(tst[cross]).median(1)
    out['max_within_abs_dpsi'] = dp[within].where(tst[within]).max(1)
    testable = replicated(tst)
    div = replicated(ds) & (out['cross_equivalent'] == 0)
    con = replicated(eq) & (out['cross_DS'] == 0)
    no_within_ds = out['within_DS'] == 0
    high = no_within_ds & (eq[hh].sum(1) >= 1) & eq[mm].all(1)
    out['status_confidence'] = np.select(
        [~testable, div & high, div, con & high, con],
        ['not testable', 'divergent high', 'divergent low', 'conserved high', 'conserved low'], 'inconclusive')
    out['status_high'] = np.select([~testable, div & high, con & high],
                                   ['not testable', 'divergent', 'conserved'], 'inconclusive')
    out['status_high_low'] = np.select([~testable, div, con], ['not testable', 'divergent', 'conserved'],
                                       'inconclusive')
    out['DS_within_species'] = ~no_within_ds
    for c in cross + within:  # per-comparison status, for the supplementary table
        out[f'status_{c}'] = st[c]
    return out.reset_index()


def main():
    variant = sys.argv[1]
    out_dir = BASE + f'events/{variant}/classification/'
    os.makedirs(out_dir, exist_ok=True)
    ann = annotate(variant)
    ann.to_csv(out_dir + 'events_annotation.tsv', sep='\t', index=False)
    print(ann['event_type'].value_counts().to_dict(), 'genes:', ann['ensembl'].nunique())
    print('splice-site correspondence:', ann['splice_site_correspondence'].value_counts(dropna=False).to_dict())

    events = pd.read_csv(BASE + f'events/{variant}/effects/event_effects.tsv', sep='\t')
    status = classify(events).merge(ann[['cluster', 'symbol', 'ensembl', 'event_type', 'splice_site_correspondence']],
                                    on='cluster', how='left')
    status.to_csv(out_dir + 'event_status_by_cell.tsv', sep='\t', index=False)
    summary = pd.concat([status.groupby('cell')[f'status_{n}'].value_counts().unstack(fill_value=0)
                         .add_prefix(f'{n}_') for n in LEVELS], axis=1)
    for n in ('divergent', 'conserved'):
        summary[f'high_low_{n}_DS_within_species'] = (status['DS_within_species'] &
                                                     (status['status_high_low'] == n)).groupby(status['cell']).sum()
    any_cell = {f'{n}_{s}': status.loc[status[f'status_{n}'] == s, 'cluster'].nunique()
                for n in LEVELS for s in ('divergent', 'conserved')}
    summary.to_csv(out_dir + 'summary_by_cell.tsv', sep='\t')
    print(summary.to_string())
    print('events in >= 1 cell type:', any_cell)


if __name__ == '__main__':
    main()
