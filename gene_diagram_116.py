#!/usr/bin/env python3
"""
Panel A of the gene figures, drawn from the Ensembl 116 data: the exons at the ends of the junctions of one AS
event in human and in mouse (transcript orientation, introns compressed), with one arc per junction in the colors
of panel B (gene_panel_116.py). Exon labels: rank in the Ensembl canonical transcript (* = not in it) and coding
length (exon widths are proportional to it); exons that share a boundary (alternative 5'/3' splice sites) are drawn
as one block with the longer exon hatched.

Usage: gene_diagram_116.py <cluster> <output prefix> --colors '#..,#..'   (one per junction, in genomic order)
"""
import argparse

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

matplotlib.rcParams['font.family'] = 'Liberation Sans'
matplotlib.rcParams['pdf.fonttype'] = 42

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
NT = 0.012          # display units per nucleotide
INTRON = 2.2        # display width of every intron
MIN_EXON = 0.6


def coord(p):
    return int(str(p).split(':')[1])


def junctions(cluster):
    idx = pd.read_csv(B + 'events/id50/AS_clusters_value.txt', sep='\t', usecols=[0]).iloc[:, 0]
    rows = idx[idx.str.endswith(':' + cluster)]
    return [tuple(r.split(':')[:3]) for r in rows]


def species_exons(pairs, sp, ends):
    """Exons (acceptor, donor, rank, coding length) with a boundary at a junction end; canonical preferred."""
    e = pairs[[f'exon_{sp}', f'rank_{sp}', f'start_{sp}', f'end_{sp}', f'len_{sp}']].drop_duplicates(f'exon_{sp}')
    e = e.assign(acc=e[f'start_{sp}'].map(coord), don=e[f'end_{sp}'].map(coord))
    minus = (e['acc'] > e['don']).iloc[0]
    hit = e[e['acc'].isin(ends) | e['don'].isin(ends)]
    out = []
    for g in sorted(ends):
        h = hit[(hit['acc'] == g) | (hit['don'] == g)]
        if h.empty:
            continue
        h = h.assign(canon=h[f'rank_{sp}'] != '*').sort_values(['canon', f'len_{sp}'], ascending=[False, True])
        r = h.iloc[0]
        out.append((r['acc'], r['don'], r[f'rank_{sp}'], int(r[f'len_{sp}'])))
    return sorted(set(out), key=lambda x: (-x[0] if minus else x[0])), minus


def layout(exons):
    """Exons sharing an acceptor or a donor form one block (alternative 5'/3' splice sites); blocks in transcript
    order. Returns [(block x0, block width, [(exon, x_acceptor, x_donor)])] and the total width."""
    blocks = []
    for ex in exons:
        for bl in blocks:
            if any(ex[0] == o[0] or ex[1] == o[1] for o in bl):
                bl.append(ex)
                break
        else:
            blocks.append([ex])
    out, x = [], 0.0
    for bl in blocks:
        main = sorted(bl, key=lambda o: o[2] == '*')[0]
        core = max(main[3] * NT, MIN_EXON)
        shared_donor = len(bl) > 1 and len({o[1] for o in bl}) == 1
        side = 0 if shared_donor else 1          # the alternative site: acceptor (0) or donor (1)
        if len(bl) > 1:
            # offset of each alternative site from the canonical one, in genomic nucleotides
            off = {o: abs(o[side] - main[side]) * NT for o in bl}
            longer = {o: abs(o[side] - main[1 - side]) > abs(main[side] - main[1 - side]) for o in bl}
            ext = max(off[o] for o in bl if longer[o]) if any(longer.values()) else 0
            w = core + ext
            items = []
            for o in bl:
                d = off[o] if longer[o] else -off[o]
                if shared_donor:
                    items.append((o, x + ext - d, x + w))
                else:
                    items.append((o, x, x + core + d))
        else:
            w = core
            items = [(bl[0], x, x + w)]
        out.append((x, w, items))
        x += w + INTRON
    return out, x - INTRON


def draw_track(ax, y, exons, juncs, colors, name):
    blocks, width = layout(exons)
    xpos = {}
    for _, _, items in blocks:
        for ex, xa, xd in items:
            xpos.setdefault(ex[0], xa)
            xpos.setdefault(ex[1], xd)
    ax.plot([0, width], [y, y], color='#555555', lw=0.8, zorder=0)
    for x0, w, items in blocks:
        items = sorted(items, key=lambda i: i[2] - i[1], reverse=True)
        for k, (ex, xa, xd) in enumerate(items):
            alt = len(items) > 1 and k == 0
            ax.add_patch(Rectangle((xa, y - 0.18), xd - xa, 0.36, facecolor='white' if alt else '#7f7f7f',
                                   edgecolor='#333333', hatch='////' if alt else None, lw=0.6, zorder=2))
        main = sorted(items, key=lambda i: i[0][2] == '*')[0][0]
        lab = [f'{main[2]} ({main[3]} nt)']
        if len(items) > 1:
            d = max(abs(i[0][0] - main[0]) or abs(i[0][1] - main[1]) for i in items)
            lab.append(f'alternative sites {d} nt apart')
        ax.text(x0 + w / 2, y - 0.3, '\n'.join(lab), ha='center', va='top', fontsize=6.5, linespacing=1.2)
    for (s, e), c in zip(juncs, colors):
        a, b = sorted((xpos[s], xpos[e]))
        h = 0.35 + 0.08 * (b - a)
        t = [i / 50 for i in range(51)]
        ax.plot([a + (b - a) * u for u in t], [y + 0.18 + h * 4 * u * (1 - u) for u in t], color=c, lw=2, zorder=1)
    ax.text(-0.4, y, name, ha='right', va='center', fontsize=9, fontweight='bold')
    return width


def plot(cluster, out, colors, title):
    pts = pd.read_csv(B + 'orthology/unique_points_116_id50.txt', sep='\t', index_col=0)
    h2m = dict(zip(pts['position_h'], pts['position_m']))
    jh = junctions(cluster)
    jm = [(coord(h2m[f'{c}:{s}']), coord(h2m[f'{c}:{e}'])) for c, s, e in jh]
    jh = [(int(s), int(e)) for _, s, e in jh]
    gene = pts.loc[pts['position_h'] == f'{junctions(cluster)[0][0]}:{jh[0][0]}', 'ensembl'].iloc[0]
    pairs = pd.read_csv(B + 'orthology/exon_pairs_116.tsv', sep='\t')
    pairs = pairs[pairs['ensembl'] == gene]
    fig, ax = plt.subplots(figsize=(8.5, 2.6))
    widths = []
    for y, sp, jj in ((1.7, 'h', jh), (0.0, 'm', jm)):
        ends = {p for j in jj for p in j}
        ex, _ = species_exons(pairs, sp, ends)
        widths.append(draw_track(ax, y, ex, jj, colors, 'Human' if sp == 'h' else 'Mouse'))
    ax.set_xlim(-1.6, max(widths) + 0.5)
    ax.set_ylim(-0.95, 3.0)
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=10)
    fig.savefig(out + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(out + '.pdf', bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cluster')
    ap.add_argument('out_prefix')
    ap.add_argument('--colors', required=True)
    ap.add_argument('--title')
    a = ap.parse_args()
    plot(a.cluster, a.out_prefix, a.colors.split(','), a.title)


if __name__ == '__main__':
    main()
