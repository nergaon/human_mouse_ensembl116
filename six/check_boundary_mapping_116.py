#!/usr/bin/env python3
"""
Boundary-mapping check of the differentially spliced (DS) events (TMSB4X-type artifact).

A junction of an event is defined by two exon boundaries in each species (human: the human junction; mouse: the
orthologous mouse boundaries, Table S7). If one boundary was mapped to a minor exon variant that the cells do not use,
the junction has (almost) no reads in that species, while the reads go to a 'near-miss' junction: the same other end
and an end a few nucleotides away (the boundary actually used). This makes an apparent human-mouse difference.

For every junction of every DS event and each species (immune samples of the six datasets):
  reads      = reads of the junction (AS_clusters_value.txt, the orthologous counts used in the analysis)
  near-miss  = the junction with the most reads in the raw species junction table (all samples of the immune studies)
               that shares one end with the junction and has the other end 1-30 nt away, and is not a junction of
               the event
Flag: near-miss reads >= 100 and reads < 5% of near-miss reads, and the junction is used in the other species
(>= 10% of the event reads there).
Output: events/id50_6ds_anova/revision/boundary_mapping_check.tsv (one row per flagged junction) and a summary.
"""
import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
ED = B + 'ensembl116/events/id50_6ds_anova/'
RAW = {
    'human': [B + f'{s}/leafcutter_0.2.9/{s}_JSR_junction_counts_no_strand.tsv'
              for s in ('GSE60424', 'GSE115736', 'GSE118165')],
    'mouse': [B + f'{s}/leafcutter_0.2.9/{s}_JSR_junction_counts_no_strand.tsv'
              for s in ('GSE287527', 'GSE313371', 'GSE261448', 'GSE184496', 'PRJNA797236')] +
             [B + f'ensembl116/m3/{s}/leafcutter_0.2.9/{s}_JSR_junction_counts_no_strand.tsv'
              for s in ('GSE215775', 'GSE307667', 'GSE263132', 'GSE236819', 'GSE254202', 'GSE231772', 'GSE134322',
                        'GSE284441', 'GSE180053')],
}
MAX_DIST, MIN_NEAR, MAX_FRAC, MIN_OTHER = 30, 100, 0.05, 0.10


def raw_near(species, ends):
    """Summed raw counts of the junctions with an end within MAX_DIST of one of ends {(chrom, pos)}."""
    keys = {(c, p + d) for c, p in ends for d in range(-MAX_DIST, MAX_DIST + 1)}
    tot = {}
    for f in RAW[species]:
        for ch in pd.read_csv(f, sep='\t', index_col=0, chunksize=200000):
            p = ch.index.str.split(':', expand=True)
            c, s, e = p.get_level_values(0), p.get_level_values(1).astype(int), p.get_level_values(2).astype(int)
            keep = [(a, b) in keys or (a, x) in keys for a, b, x in zip(c, s, e)]
            sub = ch[keep].sum(axis=1)
            for j, v in zip([f'{a}:{b}:{x}' for a, b, x in zip(c[keep], s[keep], e[keep])], sub.values):
                tot[j] = tot.get(j, 0) + v
    return pd.Series(tot, dtype=float)


def main():
    st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
    ds = st.loc[st['status_high'] == 'divergent']
    ds_cells = ds.groupby('cluster')['cell'].apply(lambda x: ','.join(sorted(x)))
    gene = st.drop_duplicates('cluster').set_index('cluster')['symbol']
    cnt = pd.read_csv(ED + 'AS_clusters_value.txt', sep='\t', index_col=0)
    clu = cnt.index.str.rsplit(':', n=1).str[1]
    cnt = cnt[clu.isin(ds_cells.index)]
    clu = cnt.index.str.rsplit(':', n=1).str[1]
    des = pd.read_csv(ED + 'sample_design.tsv', sep='\t', index_col=0)
    des = des[(des['cell'] != 'Fibro') & des.index.isin(cnt.columns)]
    reads = pd.DataFrame({sp: cnt[des.index[des['dataset'].str[0] == sp[0].upper()]].sum(axis=1)
                          for sp in ('human', 'mouse')})
    share = reads / reads.groupby(clu).transform('sum')
    # junction coordinates in each species
    pts = pd.read_csv(B + 'ensembl116/orthology/unique_points_116_id50_6ds.txt', sep='\t', index_col=0)
    h2m = dict(zip(pts['position_h'], pts['position_m']))
    coords = {}
    for j in cnt.index:
        c, s, e, _ = j.split(':')
        ms, me = h2m[f'{c}:{s}'], h2m[f'{c}:{e}']
        mc = ms.split(':')[0]
        a, b = sorted([int(ms.split(':')[1]), int(me.split(':')[1])])
        coords[j] = {'human': (c, int(s), int(e)), 'mouse': (mc, a, b)}
    rows = []
    for sp in ('human', 'mouse'):
        other = 'mouse' if sp == 'human' else 'human'
        ends = {(c, p) for v in coords.values() for c, p in ((v[sp][0], v[sp][1]), (v[sp][0], v[sp][2]))}
        raw = raw_near(sp, ends)
        rc = raw.index.str.split(':', expand=True)
        rdf = pd.DataFrame({'chrom': rc.get_level_values(0), 's': rc.get_level_values(1).astype(int),
                            'e': rc.get_level_values(2).astype(int), 'reads': raw.values}, index=raw.index)
        for cl, js in pd.Series(cnt.index, index=clu).groupby(level=0):
            own = {coords[j][sp] for j in js}
            for j in js:
                c, s, e = coords[j][sp]
                near = rdf[(rdf['chrom'] == c) & (((rdf['s'] == s) & (rdf['e'] != e) & ((rdf['e'] - e).abs() <= MAX_DIST)) |
                                                  ((rdf['e'] == e) & (rdf['s'] != s) & ((rdf['s'] - s).abs() <= MAX_DIST)))]
                near = near[[(c, a, b) not in own for a, b in zip(near['s'], near['e'])]]
                if near.empty:
                    continue
                n = near.sort_values('reads').iloc[-1]
                r = reads.loc[j, sp]
                if n['reads'] >= MIN_NEAR and r < MAX_FRAC * n['reads'] and share.loc[j, other] >= MIN_OTHER:
                    rows.append(dict(cluster=cl, gene=gene.get(cl), DS_cells=ds_cells[cl], species=sp,
                                     junction_human=j.rsplit(':', 1)[0], junction_in_species=f'{c}:{s}:{e}',
                                     reads=int(r), near_miss=n.name, near_miss_reads=int(n['reads']),
                                     distance_nt=int(abs(n['s'] - s) + abs(n['e'] - e)),
                                     share_other_species=round(share.loc[j, other], 3)))
    out = pd.DataFrame(rows)
    out.to_csv(ED + 'revision/boundary_mapping_check.tsv', sep='\t', index=False)
    n_ds = ds_cells.size
    ev = out['cluster'].nunique() if len(out) else 0
    print(f'DS events: {n_ds}; flagged: {ev} ({100 * ev / n_ds:.1f}%)')
    if len(out):
        print(out.groupby('species')['cluster'].nunique().to_string())
        print(out.drop_duplicates('cluster')['DS_cells'].str.count(',').add(1).value_counts().sort_index().to_string())
        print(out.sort_values('near_miss_reads', ascending=False).head(40).to_string(index=False))


if __name__ == '__main__':
    main()
