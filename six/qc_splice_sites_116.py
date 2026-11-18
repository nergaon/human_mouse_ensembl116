#!/usr/bin/env python3
"""
Splice-site quality control of the alternative splicing events of one orthology variant (run after build_events_6ds.py
and before LeafCutter / the classification).

An apparent human-mouse difference can arise when the two species use splice sites a few nucleotides apart and the
event contains the site of one species and the unused, orthologous position in the other. For every junction of every
event and each species (immune samples of the six datasets, sample_design.tsv):
  reads      = reads of the junction (AS_clusters_value.txt)
  near-miss  = the junction with the most reads in the raw junction table of the species (immune samples of the design)
               that shares one end with the junction, has the other end 1-30 nt away and is not part of the event
Flag: near-miss reads >= 100, junction reads < 5% of the near-miss reads, and the junction has >= 10% of the event
reads in the other species. Events with a flagged junction are removed.
Boundary support: an event is also removed if one of its boundary pairs is supported only by exon pairs compared on
fewer than 20 nt (e.g. 6 coding nt) and the liftOver projection of the human boundary does not land on the mouse
partner (<= 2 nt) -> qc_boundary_support.tsv.
Output: events/<variant>_6ds/qc_splice_sites.tsv (flagged junctions), AS_clusters_value_preqc.txt (all events) and
        AS_clusters_value.txt (events that pass).
Usage: qc_splice_sites_116.py <variant>
"""
import os
import sys

import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
E = B + 'ensembl116/'
RAW = {
    'human': [B + f'{s}/leafcutter_0.2.9/{s}_JSR_junction_counts_no_strand.tsv'
              for s in ('GSE60424', 'GSE115736', 'GSE118165')],
    'mouse': [B + f'{s}/leafcutter_0.2.9/{s}_JSR_junction_counts_no_strand.tsv'
              for s in ('GSE287527', 'GSE313371', 'GSE261448', 'GSE184496', 'PRJNA797236')] +
             [E + f'm3/{s}/leafcutter_0.2.9/{s}_JSR_junction_counts_no_strand.tsv'
              for s in ('GSE215775', 'GSE307667', 'GSE263132', 'GSE236819', 'GSE254202', 'GSE231772', 'GSE134322',
                        'GSE284441', 'GSE180053')],
}
MAX_DIST, MIN_NEAR, MAX_FRAC, MIN_OTHER = 30, 100, 0.05, 0.10
MIN_CMP = 20   # boundary support: compared length (nt) of the best supporting exon pair


def raw_near(species, ends, samples):
    """Summed raw counts (design samples) of the junctions with an end within MAX_DIST of one of ends."""
    keys = {(c, p + d) for c, p in ends for d in range(-MAX_DIST, MAX_DIST + 1)}
    tot = {}
    for f in RAW[species]:
        cols = pd.read_csv(f, sep='\t', nrows=0).columns
        study = os.path.basename(f).split('_JSR_')[0]
        use = [c for c in cols[1:] if f'{c}_{study}' in samples]   # design sample = raw column + '_' + study
        if not use:
            print('no design samples in', os.path.basename(f), '- using all samples', flush=True)
            use = list(cols[1:])
        for ch in pd.read_csv(f, sep='\t', index_col=0, usecols=[cols[0]] + use, chunksize=200000):
            p = ch.index.str.split(':', expand=True)
            c, s, e = p.get_level_values(0), p.get_level_values(1).astype(int), p.get_level_values(2).astype(int)
            keep = [(a, b) in keys or (a, x) in keys for a, b, x in zip(c, s, e)]
            sub = ch[keep].sum(axis=1)
            for j, v in zip([f'{a}:{b}:{x}' for a, b, x in zip(c[keep], s[keep], e[keep])], sub.values):
                tot[j] = tot.get(j, 0) + v
    return pd.Series(tot, dtype=float)


def main():
    v = sys.argv[1]
    d = E + f'events/{v}_6ds/'
    src = d + 'AS_clusters_value_preqc.txt'
    if not os.path.exists(src):
        os.rename(d + 'AS_clusters_value.txt', src)
    cnt = pd.read_csv(src, sep='\t', index_col=0)
    clu = cnt.index.str.rsplit(':', n=1).str[1]
    des = pd.read_csv(d + 'sample_design.tsv', sep='\t', index_col=0)
    des = des[(des['cell'] != 'Fibro') & des.index.isin(cnt.columns)]
    reads = pd.DataFrame({sp: cnt[des.index[des['dataset'].str[0] == sp[0].upper()]].sum(axis=1)
                          for sp in ('human', 'mouse')})
    share = reads / reads.groupby(clu).transform('sum')
    pts = pd.read_csv(E + f'orthology/unique_points_116_{v}.txt', sep='\t', index_col=0)
    h2m = dict(zip(pts['position_h'], pts['position_m']))
    coords = {}
    for j in cnt.index:
        c, s, e, _ = j.split(':')
        ms, me = h2m[f'{c}:{s}'], h2m[f'{c}:{e}']
        a, b = sorted([int(ms.split(':')[1]), int(me.split(':')[1])])
        coords[j] = {'human': (c, int(s), int(e)), 'mouse': (ms.split(':')[0], a, b)}
    rows = []
    for sp in ('human', 'mouse'):
        other = 'mouse' if sp == 'human' else 'human'
        ends = {(x[sp][0], p) for x in coords.values() for p in x[sp][1:]}
        raw = raw_near(sp, ends, set(des.index))
        rc = raw.index.str.split(':', expand=True)
        rdf = pd.DataFrame({'chrom': rc.get_level_values(0), 's': rc.get_level_values(1).astype(int),
                            'e': rc.get_level_values(2).astype(int), 'reads': raw.values}, index=raw.index)
        by_chrom = {c: g for c, g in rdf.groupby('chrom')}
        for cl, js in pd.Series(cnt.index, index=clu).groupby(level=0):
            own = {coords[j][sp] for j in js}
            for j in js:
                c, s, e = coords[j][sp]
                g = by_chrom.get(c)
                if g is None:
                    continue
                near = g[((g['s'] == s) & (g['e'] != e) & ((g['e'] - e).abs() <= MAX_DIST)) |
                         ((g['e'] == e) & (g['s'] != s) & ((g['s'] - s).abs() <= MAX_DIST))]
                near = near[[(c, a, b) not in own for a, b in zip(near['s'], near['e'])]]
                if near.empty:
                    continue
                n = near.sort_values('reads').iloc[-1]
                r = reads.loc[j, sp]
                if n['reads'] >= MIN_NEAR and r < MAX_FRAC * n['reads'] and share.loc[j, other] >= MIN_OTHER:
                    rows.append(dict(cluster=cl, species=sp, junction_human=j.rsplit(':', 1)[0],
                                     junction_in_species=f'{c}:{s}:{e}', reads=int(r), near_miss=n.name,
                                     near_miss_reads=int(n['reads']),
                                     distance_nt=int(abs(n['s'] - s) + abs(n['e'] - e)),
                                     share_other_species=round(share.loc[j, other], 3)))
    out = pd.DataFrame(rows, columns=['cluster', 'species', 'junction_human', 'junction_in_species', 'reads',
                                      'near_miss', 'near_miss_reads', 'distance_nt', 'share_other_species'])
    out.to_csv(d + 'qc_splice_sites.tsv', sep='\t', index=False)
    # 2. boundary support: a boundary pair of the event must be supported by an exon pair compared on at least
    # MIN_CMP nt, or by the liftOver projection of the human boundary landing on the mouse partner (<= 2 nt)
    unsup = unsupported_boundaries(v, pts, set(f'{c}:{p}' for j in cnt.index for c, p in
                                              [(j.split(':')[0], j.split(':')[1]), (j.split(':')[0], j.split(':')[2])]))
    unsup.to_csv(d + 'qc_boundary_support.tsv', sep='\t', index=False)
    bad_b = set(unsup['position_h'])
    hit = pd.Series([any(f"{j.split(':')[0]}:{q}" in bad_b for q in j.split(':')[1:3]) for j in cnt.index],
                    index=cnt.index)
    bad_support = set(clu[hit])
    bad_near = set(out['cluster'])
    bad = bad_near | bad_support
    keep = ~clu.isin(bad)
    cnt[keep].to_csv(d + 'AS_clusters_value.txt', sep='\t')
    print(f'{v}: events {clu.nunique()}, removed {len(bad)} ({100 * len(bad) / clu.nunique():.1f}%): '
          f'near-miss junction {len(bad_near)} (human {out.loc[out.species == "human", "cluster"].nunique()}, '
          f'mouse {out.loc[out.species == "mouse", "cluster"].nunique()}), unsupported boundary {len(bad_support)} '
          f'({len(bad_support - bad_near)} not already removed; {len(bad_b)} boundaries); kept {clu[keep].nunique()} '
          f'events, {keep.sum()} junctions')


def unsupported_boundaries(v, pts, used):
    """Boundary pairs (used by the events) supported only by exon pairs compared on < MIN_CMP nt (len_h/len_m of
    exon_pairs_116.tsv: coding length for coding exons) whose human boundary does not lift onto the mouse partner."""
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'liftover'))
    from orthology_points_116 import VARIANTS
    from liftover_boundaries import Chain, CHAIN
    min_id, rbh = VARIANTS[v]
    p = pd.read_csv(E + 'orthology/exon_pairs_116.tsv', sep='\t')
    p = p[p[['%id_human', '%id_mouse']].min(axis=1) >= min_id]
    if rbh:
        p = p[p['reciprocal_best']]
    n = p.groupby('ensembl').size()
    p = p[p['ensembl'].isin(n[n > 1].index)]
    p['cmp'] = p[['len_h', 'len_m']].min(axis=1)
    best = pd.concat([p[['start_h', 'start_m', 'cmp']].set_axis(['h', 'm', 'cmp'], axis=1),
                      p[['end_h', 'end_m', 'cmp']].set_axis(['h', 'm', 'cmp'], axis=1)]).groupby(['h', 'm'])['cmp'].max()
    q = pts[pts['position_h'].isin(used)].copy()
    q['max_compared_nt'] = [best.get((h, m)) for h, m in zip(q['position_h'], q['position_m'])]
    short = q[q['max_compared_nt'].fillna(0) < MIN_CMP].copy()   # NaN: partner set by the liftOver rule
    chain = Chain(CHAIN, set(short['position_h'].str.split(':').str[0]))
    lifted, ok = [], []
    for h, m in zip(short['position_h'], short['position_m']):
        c, pos = h.split(':')
        lf = chain.lift(c, int(pos))
        lifted.append(None if lf is None else f'{lf[0]}:{lf[1]}')
        ok.append(lf is not None and lf[0] == m.split(':')[0] and abs(lf[1] - int(m.split(':')[1])) <= 2)
    short['lifted'], short['lift_on_partner'] = lifted, ok
    return short.loc[~short['lift_on_partner'], ['position_h', 'position_m', 'symbol', 'max_compared_nt', 'lifted']]


if __name__ == '__main__':
    main()
