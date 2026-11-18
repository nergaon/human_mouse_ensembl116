#!/usr/bin/env python3
"""
LiftOver check of all human-mouse orthologous exon boundary pairs (id50, six datasets).

Ambiguous boundary pair: the human boundary H or its mouse partner M (orthology/unique_points_116_id50_6ds.txt) has
another annotated splice site of the same type (donor / acceptor, Ensembl 116, all transcripts; transcript starts
and ends are not splice sites) within 30 nt in the same gene region.

For every ambiguous pair, H is lifted to mm39 (UCSC hg38ToMm39 chain, best-scoring chain that contains the base) and
the mouse candidates (annotated mouse sites of the same type within 30 nt of M, M included) are compared with the
lifted position L:
  agree           M is the candidate nearest to L (or tied with it)
  disagree        another candidate is nearer to L; it is proposed as the partner if it is within MAX_SNAP nt of L
  no_lift / far   H cannot be lifted, or L is more than 30 nt from every candidate (alignment does not help)
In the same way, for a human site with a near annotated human site, each human candidate is lifted and the one that
lands on M is reported (human-side ambiguity, e.g. NAGNAG acceptors).

Decision rule (column 'decision'): switch the mouse partner only when the liftover points to another candidate within
MAX_SNAP nt that is a site of the Ensembl canonical transcript while the current partner is not. Liftover alone is not
used, because chain alignments place gaps ambiguously near splice sites and moved many canonical partners to minor
sites (evaluate_changes.py). A switch to a mouse site that is already paired with another human boundary is not
made (one-to-one pairing). The rule uses no RNA-seq data; evaluate_changes.py checks it with the reads.
Outputs (ensembl116/liftover_check/):
  ambiguous_boundaries.tsv               every ambiguous pair with its liftover decision
  unique_points_116_id50_6ds_liftover.txt  the boundary pairs with the proposed partners replaced
  flagged_events_liftover.tsv            the 12 events of the read-based check (boundary_mapping_check.tsv) and what
                                         the liftover decides for their boundaries
"""
import bisect
import gzip
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
OUT = B + 'liftover_check/'
CHAIN = OUT + 'ref/hg38ToMm39.over.chain.gz'
GTF = {'h': B + 'ref/Homo_sapiens.GRCh38.116.gtf.gz', 'm': B + 'ref/Mus_musculus.GRCm39.116.gtf.gz'}
# usage: liftover_boundaries.py [variant] [input points]  (default id50, orthology/unique_points_116_<variant>.txt;
# after the corrected points are installed, the input is the backup in pre_liftover_20261118/orthology/)
VARIANT = sys.argv[1] if len(sys.argv) > 1 else 'id50'
POINTS = sys.argv[2] if len(sys.argv) > 2 else B + f'orthology/unique_points_116_{VARIANT}.txt'
EVENTS = B + f'events/{VARIANT}_6ds/AS_clusters_value.txt'
if VARIANT != 'id50':
    OUT = OUT + VARIANT + '/'
FLAGGED = B + 'events/id50_6ds_anova/revision/boundary_mapping_check.tsv'
WINDOW, MAX_SNAP = 30, 2


def splice_sites(path):
    """{(chrom, pos): (type, strand, canonical)}: donors (3' exon end) and acceptors (5' exon end) of internal exon
    ends, 1-based exon boundary base, as in unique_points."""
    ex = defaultdict(list)
    canon = set()
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split('\t', 9)
            if p[2] != 'exon':
                continue
            t = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
            c = 'chrM' if p[0] == 'MT' else p[0] if p[0].startswith('chr') else 'chr' + p[0]
            ex[t].append((c, int(p[3]), int(p[4]), p[6]))
            if 'Ensembl_canonical' in p[8]:
                canon.add(t)
    sites = {}
    for t, es in ex.items():
        es.sort(key=lambda e: e[1], reverse=es[0][3] == '-')   # transcript order
        for i, (c, s, e, st) in enumerate(es):
            five, three = (s, e) if st == '+' else (e, s)
            if i > 0:
                k = (c, five)
                old = sites.get(k)
                sites[k] = ('acceptor', st, (old is not None and old[2]) or t in canon)
            if i < len(es) - 1:
                k = (c, three)
                old = sites.get(k)
                sites[k] = ('donor', st, (old is not None and old[2]) or t in canon)
    return sites


def index_sites(sites):
    by = defaultdict(list)
    for (c, p), (ty, st, _) in sites.items():
        by[(c, st, ty)].append(p)
    return {k: np.array(sorted(v)) for k, v in by.items()}


def near(idx, c, st, ty, p, w=WINDOW):
    a = idx.get((c, st, ty))
    if a is None:
        return []
    i, j = np.searchsorted(a, p - w), np.searchsorted(a, p + w, side='right')
    return list(a[i:j])


class Chain:
    """Point liftOver: for a base, the best-scoring chain whose aligned blocks contain it."""

    def __init__(self, path, chroms):
        self.blocks = defaultdict(list)   # tName -> [(tStart, tEnd, score, qName, qStrand, qSize, qStart)]
        with gzip.open(path, 'rt') as f:
            keep = False
            for line in f:
                if line.startswith('chain'):
                    h = line.split()
                    score, tn, ts, qn, qsize, qst, qs = float(h[1]), h[2], int(h[5]), h[7], int(h[8]), h[9], int(h[10])
                    keep = tn in chroms
                    t, q = ts, qs
                    continue
                if not keep or not line.strip():
                    continue
                v = line.split()
                size = int(v[0])
                self.blocks[tn].append((t, t + size, score, qn, qst, qsize, q))
                if len(v) == 3:
                    t += size + int(v[1])
                    q += size + int(v[2])
        self.starts = {}
        for tn in self.blocks:
            self.blocks[tn].sort()
            self.starts[tn] = [b[0] for b in self.blocks[tn]]
        self.maxlen = {tn: max(b[1] - b[0] for b in bl) for tn, bl in self.blocks.items()}

    def lift(self, chrom, pos1):
        """1-based position -> (qChrom, 1-based position, qStrand) or None."""
        if chrom not in self.blocks:
            return None
        p0 = pos1 - 1
        bl, st = self.blocks[chrom], self.starts[chrom]
        i = bisect.bisect_right(st, p0) - 1
        best = None
        while i >= 0 and bl[i][0] > p0 - self.maxlen[chrom] - 1:
            ts, te, score, qn, qst, qsize, qs = bl[i]
            if ts <= p0 < te and (best is None or score > best[0]):
                q0 = qs + (p0 - ts)
                best = (score, qn, (qsize - q0) if qst == '-' else q0 + 1, qst)
            i -= 1
        return None if best is None else best[1:]


def main():
    pts = pd.read_csv(POINTS, sep='\t', index_col=0)
    ids = pd.read_csv(EVENTS, sep='\t', usecols=[0]).iloc[:, 0]
    parts = ids.str.split(':', expand=True)
    used = set(parts[0] + ':' + parts[1]) | set(parts[0] + ':' + parts[2])
    ev = pts.copy()   # all orthologous boundary pairs (the events are rebuilt from the corrected pairs)
    print('boundary pairs:', len(ev), '| used by the current events:', ev['position_h'].isin(used).sum())
    sh, sm = splice_sites(GTF['h']), splice_sites(GTF['m'])
    ih, im = index_sites(sh), index_sites(sm)
    hc = ev['position_h'].str.split(':', expand=True)
    mc = ev['position_m'].str.split(':', expand=True)
    ev['hc'], ev['hp'], ev['mc'], ev['mp'] = hc[0], hc[1].astype(int), mc[0], mc[1].astype(int)
    rows = []
    for k, r in ev.iterrows():
        h, m = sh.get((r['hc'], r['hp'])), sm.get((r['mc'], r['mp']))
        if h is None or m is None:
            rows.append(dict(point=k, status='not_annotated_site'))
            continue
        cand_m = near(im, r['mc'], m[1], m[0], r['mp'])
        cand_h = near(ih, r['hc'], h[1], h[0], r['hp'])
        if len(cand_m) < 2 and len(cand_h) < 2:
            continue
        rows.append(dict(point=k, type=h[0], strand_h=h[1], strand_m=m[1], cand_m=cand_m, cand_h=cand_h))
    amb = pd.DataFrame(rows).set_index('point')
    amb = amb[amb.get('status').isna()] if 'status' in amb else amb
    print('ambiguous boundary pairs:', len(amb))
    chain = Chain(CHAIN, set(ev['hc']))
    out = []
    for k, a in amb.iterrows():
        r = ev.loc[k]
        L = chain.lift(r['hc'], r['hp'])
        rec = dict(point=k, symbol=r['symbol'], position_h=r['position_h'], position_m=r['position_m'],
                   type=a['type'], n_mouse_candidates=len(a['cand_m']), n_human_candidates=len(a['cand_h']),
                   mouse_canonical=sm[(r['mc'], r['mp'])][2], human_canonical=sh[(r['hc'], r['hp'])][2])
        if L is None or L[0] != r['mc']:
            rec.update(lifted=None if L is None else f'{L[0]}:{L[1]}', status='no_lift')
        else:
            d = {p: abs(p - L[1]) for p in a['cand_m']}
            best = min(d.values())
            rec.update(lifted=f'{L[0]}:{L[1]}', dist_chosen=d.get(r['mp'], abs(r['mp'] - L[1])), dist_best=best)
            if best > WINDOW:
                rec['status'] = 'far'
            elif d.get(r['mp'], 1e9) <= best:
                rec['status'] = 'agree'
            else:
                alt = min(d, key=d.get)
                rec.update(status='disagree', proposed_m=f"{r['mc']}:{alt}", proposed_shift=alt - r['mp'],
                           proposed_canonical=sm[(r['mc'], alt)][2], snap_ok=best <= MAX_SNAP)
        # human side: which human candidate lifts onto the chosen mouse site
        if len(a['cand_h']) > 1:
            hits = [p for p in a['cand_h'] if (x := chain.lift(r['hc'], p)) and x[0] == r['mc'] and abs(x[1] - r['mp']) <= MAX_SNAP]
            rec['human_candidates_lifting_to_M'] = ','.join(map(str, hits))
            rec['human_side'] = ('agree' if r['hp'] in hits else 'disagree' if hits else 'none_lands')
        out.append(rec)
    res = pd.DataFrame(out)
    res['in_current_events'] = res['position_h'].isin(used)
    os.makedirs(OUT, exist_ok=True)
    res.to_csv(OUT + 'ambiguous_boundaries.tsv', sep='\t', index=False)
    print(res['status'].value_counts().to_string())
    if 'human_side' in res:
        print('human side:', res['human_side'].value_counts().to_dict())
    dis_all = res[(res['status'] == 'disagree') & (res['snap_ok'] == True)]
    print('liftover disagrees (mouse side, lift within', MAX_SNAP, 'nt):', len(dis_all))
    # decision rule: switch only when the liftover partner is the canonical site and the current one is not
    # (liftover alone moved many correct canonical partners to minor sites: changed_junction_reads.tsv, first run)
    dis = dis_all[(dis_all['proposed_canonical'] == True) & (dis_all['mouse_canonical'] == False)]
    # keep the pairing one-to-one: no switch to a mouse site that is already the partner of another human boundary
    taken = set(pts['position_m'])
    print('switches to a mouse site already paired (not done):', dis['proposed_m'].isin(taken).sum())
    dis = dis[~dis['proposed_m'].isin(taken)]
    res['decision'] = 'keep'
    res.loc[dis.index, 'decision'] = 'switch'
    res.to_csv(OUT + 'ambiguous_boundaries.tsv', sep='\t', index=False)
    print('partner changes (liftover partner canonical, current not):', len(dis))
    # proposed mapping
    new = pts.copy()
    new.loc[dis['point'], 'position_m'] = dis['proposed_m'].values
    os.makedirs(OUT, exist_ok=True)
    new.to_csv(OUT + f'unique_points_116_{VARIANT}_liftover.txt', sep='\t')
    # events affected
    jp = pd.DataFrame({'clu': ids.str.rsplit(':', n=1).str[1], 's': parts[0] + ':' + parts[1],
                       'e': parts[0] + ':' + parts[2]})
    chpos = set(dis['position_h'])
    print('events with a changed mouse partner:', jp.loc[jp['s'].isin(chpos) | jp['e'].isin(chpos), 'clu'].nunique())
    if VARIANT != 'id50':
        return
    # the 12 events of the read-based check
    fl = pd.read_csv(FLAGGED, sep='\t')
    rows = []
    for _, f in fl.iterrows():
        c, s, e = f['junction_human'].split(':')
        for p in (f'{c}:{s}', f'{c}:{e}'):
            x = res[res['position_h'] == p]
            rows.append(dict(cluster=f['cluster'], gene=f['gene'], species_without_reads=f['species'],
                             near_miss=f['near_miss'], position_h=p,
                             position_m=pts.loc[pts['position_h'] == p, 'position_m'].iloc[0],
                             liftover=x['status'].iloc[0] if len(x) else 'not ambiguous',
                             lifted=x['lifted'].iloc[0] if len(x) else None,
                             proposed_m=x['proposed_m'].iloc[0] if len(x) and 'proposed_m' in x else None,
                             human_side=x['human_side'].iloc[0] if len(x) and 'human_side' in x else None,
                             human_candidates_lifting_to_M=x['human_candidates_lifting_to_M'].iloc[0]
                             if len(x) and 'human_candidates_lifting_to_M' in x else None))
    fr = pd.DataFrame(rows)
    fr.to_csv(OUT + 'flagged_events_liftover.tsv', sep='\t', index=False)
    print(fr.to_string(index=False))


if __name__ == '__main__':
    main()
