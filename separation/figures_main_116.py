#!/usr/bin/env python3
"""
Main data figures of the revision (Ensembl 116, id50).

  Figure 2   classification: A human-human, B mouse-mouse, C human-mouse (+ fibroblasts) status
             of testable events (mean over dataset pairs); D divergent / conserved events of high and low
             confidence per cell type; E/F UpSet of high-confidence divergent / conserved events; G heatmap of
             events with a high-confidence class in >= 3 immune cell types (+ fibroblasts)
  Figure S2  A-C as Figure 2 as % of all orthologous events; E/F with high- and low-confidence events;
             G events with a high-confidence class in exactly 2 immune cell types
  Figure 3   effect sizes: A cumulative distribution of event max|dPSI| on events testable in every
             comparison, one line per dataset pair; B Pearson r of junction PSI per dataset pair;
             C PSI scatter for one human-human, mouse-mouse and human-mouse pair (NK);
             D high / high+low confidence fractions by minor-isoform usage
  Figure 6   fibroblasts: A-D fibroblast status of immune species-divergent / conserved events (high or low
             confidence) and immune class of fibroblast DS / equivalent events (per cell type)
Output: ensembl116/revision_separation/figures/Figure_<n>.pdf / .tif
"""
import itertools
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from scipy import stats

matplotlib.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 7, 'pdf.fonttype': 42,
                            'axes.spines.top': False, 'axes.spines.right': False})
ED = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/id50_sep/'
OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_separation/figures/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_LABEL = {'CD4T': 'T4', 'CD8T': 'T8', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo', 'Fibro': 'Fb'}
STATUS_COLORS = {'DS': '#74c476', 'equivalent': '#9ecae1', 'inconclusive': '#d9d9d9'}
STATUS_LABEL = {'DS': 'Differentially spliced', 'equivalent': 'Equivalent (±0.1)', 'inconclusive': 'Inconclusive'}
ALL_COLORS = {**STATUS_COLORS, 'not_testable': '#f7f7f7'}   # Figure S2 A-C: all events, incl. not testable
ALL_LABEL = {**STATUS_LABEL, 'not_testable': 'Not testable'}
STATUS = 'status_m0.1'  # per-comparison status: DS, equivalent within +-0.1, inconclusive
HC_COLORS = {'divergent': '#238b45', 'conserved': '#2171b5', 'inconclusive': '#d9d9d9'}
CONF_COLORS = {'divergent high': '#238b45', 'divergent low': '#a1d99b', 'conserved high': '#2171b5',
               'conserved low': '#9ecae1', 'inconclusive': '#d9d9d9'}
NOT_TESTABLE = '#f7f7f7'
HIGH_HATCH = '....'
HATCH_COLOR = '#333333'
matplotlib.rcParams['hatch.linewidth'] = 0.6
KIND_COLORS = {'HH': '#2a78d6', 'cross': '#eb6834', 'MM': '#1baf7a'}
KIND_LABEL = {'HH': 'Human-Human', 'cross': 'Human-Mouse', 'MM': 'Mouse-Mouse'}


def kind(c):
    if c.startswith('matched_'):
        return 'matched'
    if c == 'FH_FM':
        return 'fibro'
    a, b = c.split('_')
    return 'cross' if a[0] != b[0] else ('HH' if a[0] == 'H' else 'MM')


def save(fig, name):
    fig.savefig(OUT + name + '.pdf', bbox_inches='tight')
    fig.savefig(OUT + name + '.tif', bbox_inches='tight', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
    print('wrote', name)


def label(ax, letter, x=-0.18, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')


# ------------------------------------------------------------------ Figure 2 / S2
def status_bars(ax, ev, k, title, fibro=None, all_events=False):
    """Median over the dataset pairs (as in Table 1) of the % of testable events (Figure 2) or, with all_events, of all
    orthologous events, including those not testable (Figure S2)."""
    colors = ALL_COLORS if all_events else STATUS_COLORS
    rows = []
    groups = [(CELL_LABEL[c], ev[(ev['kind'] == k) & (ev['cell'] == c)]) for c in CELLS] if k else []
    if fibro is not None:
        groups.append(('Fb', fibro))
    for lab, g in groups:
        if not all_events:
            g = g[g[STATUS] != 'not_testable']
        frac = g.groupby('comparison')[STATUS].value_counts(normalize=True).unstack(fill_value=0).median()
        rows.append(frac.rename(lab))
    df = pd.DataFrame(rows).reindex(columns=list(colors)).fillna(0) * 100
    counts = pd.DataFrame([(g if all_events else g[g[STATUS] != 'not_testable']).groupby('comparison')[STATUS]
                           .value_counts().unstack(fill_value=0).median().rename(lab) for lab, g in groups]) \
        .reindex(columns=list(colors)).fillna(0)
    bottom = np.zeros(len(df))
    for s, c in colors.items():
        ax.bar(df.index, df[s], bottom=bottom, color=c, width=0.75, edgecolor='white' if s != 'not_testable' else
               '#bdbdbd', linewidth=0.5)
        for x, (b, v, k) in enumerate(zip(bottom, df[s], counts[s])):
            if v >= 7:
                ax.text(x, b + v / 2, f'{k:.0f}', ha='center', va='center', fontsize=4.8, color='#333333')
        bottom += df[s].values
    ax.set_ylim(0, 100)
    if len(df) == 1:
        ax.set_xlim(-1.5, 1.5)
    ax.set_title(title, fontsize=8)
    ax.set_ylabel('% of all orthologous events' if all_events else '% of testable events')
    ax.tick_params(axis='x', rotation=0)


def hc_bars(ax, status):
    """% of testable events per class; high confidence dotted, low confidence plain (colors as in A-C)."""
    t = status[status['status_confidence'] != 'not testable']
    df = pd.crosstab(t['cell'], t['status_confidence'], normalize='index').reindex(CELLS) * 100
    df.index = [CELL_LABEL[c] for c in df.index]
    bottom = np.zeros(len(df))
    for s in CONF_COLORS:
        v = df.get(s, pd.Series(0, index=df.index))
        c = STATUS_COLORS[{'divergent': 'DS', 'conserved': 'equivalent'}.get(s.split()[0], 'inconclusive')]
        ax.bar(df.index, v, bottom=bottom, color=c, width=0.75, hatch=HIGH_HATCH if s.endswith('high') else None,
               edgecolor=HATCH_COLOR, linewidth=0)
        ax.bar(df.index, v, bottom=bottom, width=0.75, fill=False, edgecolor='white', linewidth=0.5)
        bottom += v.values
    ax.set_ylim(0, 100)
    ax.set_ylabel('% of testable events')
    ax.set_title('Human-mouse classification', fontsize=8)


def upset(fig, spec, status, col, cls, color, title):
    sets = {c: set(status.loc[(status['cell'] == c) & (status[col] == cls), 'cluster']) for c in CELLS}
    events = set().union(*sets.values())
    combo = pd.Series({e: tuple(c for c in CELLS if e in sets[c]) for e in events}).value_counts()
    combo = combo.head(25)
    gs = spec.subgridspec(2, 1, height_ratios=[1.6, 1], hspace=0.05)
    ax1, ax2 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    x = np.arange(len(combo))
    ax1.bar(x, combo.values, color='#8c8c8c', width=0.6)
    for xi, v in zip(x, combo.values):
        ax1.text(xi, v, str(v), ha='center', va='bottom', fontsize=5)
    ax1.set_xlim(-0.6, len(combo) - 0.4)
    ax1.set_xticks([])
    ax1.set_ylabel('Events')
    ax1.set_title(title, fontsize=8)
    for yi, c in enumerate(CELLS):
        member = [c in k for k in combo.index]
        ax2.scatter(x, [yi] * len(x), s=10, c=[color if m else '#e0e0e0' for m in member], zorder=2)
    for xi, k in enumerate(combo.index):
        ys = [CELLS.index(c) for c in k]
        ax2.plot([xi, xi], [min(ys), max(ys)], color=color, lw=1, zorder=1)
    ax2.set_yticks(range(len(CELLS)))
    ax2.set_yticklabels([CELL_LABEL[c] for c in CELLS])
    ax2.set_xlim(-0.6, len(combo) - 0.4)
    ax2.set_xticks([])
    ax2.invert_yaxis()
    for s in ax2.spines.values():
        s.set_visible(False)
    return ax1


def heatmap(ax, status, col, fib, cells=(4, 5)):
    """Events with a high-confidence class in cells[0]-cells[1] immune cell types; colors as in A-C, a dot marks
    high confidence in an immune cell type; fibroblasts by their single-pair status."""
    t = status[status['status_high'].isin(['divergent', 'conserved'])]
    n = t.groupby('cluster')['cell'].nunique()
    keep = n[n.between(*cells)].index
    sub = status[status['cluster'].isin(keep)]
    w = sub.pivot_table(index='cluster', columns='cell', values=col, aggfunc='first').reindex(columns=CELLS)
    conf = sub.pivot_table(index='cluster', columns='cell', values='status_confidence',
                           aggfunc='first').reindex(index=w.index, columns=CELLS)
    conf['Fibro'] = w.index.map(fib).fillna('not testable')
    cls = {'divergent high': 1, 'divergent low': 1, 'DS': 1, 'conserved high': 2, 'conserved low': 2,
           'equivalent': 2, 'inconclusive': 3}
    m = conf.apply(lambda s: s.map(lambda v: cls.get(v, 0)))
    hi = conf.isin(['divergent high', 'conserved high'])
    # DS (high confidence in any cell type) first, then equivalent; within each, by number of cell types
    high = sub.pivot_table(index='cluster', columns='cell', values='status_high', aggfunc='first').reindex(w.index)
    first = high.apply(lambda r: 2 if (r == 'divergent').any() else 1, axis=1)
    m = m.assign(_k=first, _n=n.reindex(m.index)).sort_values(['_k', '_n'], ascending=False).drop(columns=['_k', '_n'])
    hi = hi.reindex(m.index)
    cmap = matplotlib.colors.ListedColormap([NOT_TESTABLE, STATUS_COLORS['DS'], STATUS_COLORS['equivalent'],
                                             STATUS_COLORS['inconclusive']])
    ax.imshow(m.values, aspect='auto', cmap=cmap, vmin=0, vmax=3, interpolation='none')
    for j in range(hi.shape[1]):                       # high confidence: dotted, as in D (one patch per run)
        v = np.r_[False, hi.iloc[:, j].fillna(False).to_numpy(bool), False]
        starts, ends = np.flatnonzero(~v[:-1] & v[1:]), np.flatnonzero(v[:-1] & ~v[1:])
        for a, b in zip(starts, ends):
            ax.add_patch(Rectangle((j - 0.5, a - 0.5), 1, b - a, fill=False, hatch=HIGH_HATCH,
                                   edgecolor=HATCH_COLOR, linewidth=0))
    sym = status.drop_duplicates('cluster').set_index('cluster')['symbol']
    # long heatmaps (Figure S2G): label every k-th row so that the labels stay readable
    step = int(np.ceil(len(m) / 150))
    rows = list(range(0, len(m), step))
    ax.set_yticks(rows)
    ax.set_yticklabels(sym.reindex(m.index[rows]).fillna(''), fontsize=max(2.5, min(5, 300 / max(len(rows), 1))))
    ax.tick_params(axis='y', length=1, width=0.3)
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels([CELL_LABEL[c] for c in m.columns], rotation=45, ha='left')
    ax.xaxis.tick_top()
    ax.axvline(len(CELLS) - 0.5, color='white', lw=3)
    ax.set_xlim(-0.5, m.shape[1] - 0.5)
    for s in ax.spines.values():
        s.set_visible(False)
    return len(m)


def figure2(ev, status, definition, name, cells):
    col = f'status_{definition}'
    all_events = definition == 'high_low'      # Figure S2 A-C: % of all orthologous events
    fib = ev[ev['comparison'] == 'FH_FM']
    fib_status = fib.set_index('cluster')[STATUS]
    fig = plt.figure(figsize=(7.2, 9.8))
    gs = fig.add_gridspec(4, 4, width_ratios=[1, 1, 1.25, 1.1], height_ratios=[1, 1, 1.05, 1.05], hspace=0.6,
                          wspace=0.5)
    ax = fig.add_subplot(gs[0, 0]); status_bars(ax, ev, 'HH', 'Human-Human', all_events=all_events); label(ax, 'A')
    ax = fig.add_subplot(gs[0, 1]); status_bars(ax, ev, 'MM', 'Mouse-Mouse', all_events=all_events)
    ax.set_ylabel(''); label(ax, 'B')
    ax = fig.add_subplot(gs[0, 2]); status_bars(ax, ev, 'cross', 'Human-Mouse', fib, all_events=all_events)
    ax.set_ylabel(''); label(ax, 'C')
    ax = fig.add_subplot(gs[1, 2]); hc_bars(ax, status); label(ax, 'D')
    # one legend for all panels, to the left of D
    lax = fig.add_subplot(gs[1, 0:2]); lax.axis('off')
    entries = [(Patch(facecolor=STATUS_COLORS['DS'], linewidth=0), 'Differentially spliced'),
               (Patch(facecolor=STATUS_COLORS['equivalent'], linewidth=0), 'Equivalent'),
               (Patch(facecolor=STATUS_COLORS['inconclusive'], linewidth=0), 'Inconclusive')]
    entries += [(Patch(facecolor=NOT_TESTABLE, edgecolor='#bdbdbd', linewidth=0.5),
                 'Not testable (A-C, G)' if all_events else 'Not testable (G)'),
                (Patch(facecolor='white', hatch=HIGH_HATCH, edgecolor=HATCH_COLOR, linewidth=0), 'High confidence (D, G)'),
                (Patch(facecolor='white', edgecolor='#bdbdbd', linewidth=0.5), 'Low confidence (D, G)')]
    lax.legend(*zip(*entries), fontsize=6.5, frameon=False, loc='center', handlelength=2.2, handleheight=1.4,
               labelspacing=0.7)
    a1 = upset(fig, gs[2, 0:3], status, col, 'divergent', HC_COLORS['divergent'], 'Differentially spliced')
    label(a1, 'E', x=-0.06)
    a1 = upset(fig, gs[3, 0:3], status, col, 'conserved', HC_COLORS['conserved'], 'Equivalent')
    label(a1, 'F', x=-0.06)
    ax = fig.add_subplot(gs[:, 3])
    n = heatmap(ax, status, col, fib_status, cells)
    label(ax, 'G', x=-0.5, y=1.02)
    fig.text(0.02, 0.005, name.replace('_', ' '), fontsize=8)
    save(fig, name)
    return n


# ------------------------------------------------------------------ Figure 3
def figure3(ev, jun, usage):
    e = ev[ev['kind'].isin(['HH', 'MM', 'cross'])]
    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 5, height_ratios=[1, 0.9, 1.1], hspace=0.7, wspace=0.15)
    for i, cell in enumerate(CELLS):
        ax = fig.add_subplot(gs[0, i])
        ec = e[e['cell'] == cell]
        n_comp = ec['comparison'].nunique()
        tested = ec[ec['testable']].groupby('cluster')['comparison'].nunique()
        common = tested[tested == n_comp].index
        for comp, g in ec[ec['cluster'].isin(common)].groupby('comparison'):
            x = np.sort(g['max_abs_dpsi'].dropna())
            ax.step(x, np.arange(1, len(x) + 1) / len(x), where='post', color=KIND_COLORS[kind(comp)], lw=1,
                    alpha=0.85)
        ax.axvline(0.1, color='#6b6b6b', lw=0.6, ls=':')
        ax.set_xlim(0, 0.6)
        ax.set_xticks([0, 0.25, 0.5])
        ax.set_ylim(0, 1.01)
        ax.set_title(f'{CELL_LABEL[cell]} (n={len(common)})', fontsize=7)
        ax.set_xlabel('event max |ΔPSI|')
        if i == 0:
            ax.set_ylabel('cumulative fraction of events')
            label(ax, 'A', x=-0.45)
        else:
            ax.tick_params(labelleft=False)
    fig.legend([plt.Line2D([], [], color=c, lw=1.5) for c in KIND_COLORS.values()], KIND_LABEL.values(),
               loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.95), fontsize=7)

    row = gs[1, :].subgridspec(1, 2, width_ratios=[1.25, 1], wspace=0.45)
    ax = fig.add_subplot(row[0])
    # B: same events as A (testable in all ten comparisons of the cell type)
    common_by_cell = {}
    for cell in CELLS:
        ec = e[e['cell'] == cell]
        tested = ec[ec['testable']].groupby('cluster')['comparison'].nunique()
        common_by_cell[cell] = set(tested[tested == ec['comparison'].nunique()].index)
    keep = pd.DataFrame([(c, k) for c, ks in common_by_cell.items() for k in ks], columns=['cell', 'cluster'])
    jc = jun[jun['testable']].merge(keep, on=['cell', 'cluster'])
    corr = jc.groupby(['comparison', 'cell']).apply(
        lambda g: stats.pearsonr(g['psi1'], g['psi2'])[0]).rename('r').reset_index()
    corr = corr[corr['comparison'].map(kind).isin(KIND_COLORS)]
    offs = {'HH': -0.22, 'MM': 0.0, 'cross': 0.22}
    for _, r in corr.iterrows():
        k = kind(r['comparison'])
        ax.scatter(CELLS.index(r['cell']) + offs[k] + np.random.uniform(-0.05, 0.05), r['r'], s=12,
                   color=KIND_COLORS[k], edgecolor='white', linewidth=0.3)
    ax.set_xticks(range(len(CELLS)))
    ax.set_xticklabels([CELL_LABEL[c] for c in CELLS])
    ax.set_ylabel('Pearson r of junction PSI')
    ax.set_ylim(0.65, 1.0)
    label(ax, 'B', x=-0.12)

    ax = fig.add_subplot(row[1])
    u = usage[usage['status_high'] != 'not testable'].copy()
    u['bin'] = pd.cut(u['minor_usage'], [-0.001, 0.01, 0.05, 0.1, 0.2, 0.5],
                      labels=['<1', '1-5', '5-10', '10-20', '20-50'])
    for col, ls, lab in (('status_high', '-', 'high'), ('status_high_low', '--', 'high + low')):
        div = u.groupby('bin', observed=False)[col].apply(lambda s: 100 * (s == 'divergent').mean())
        con = u.groupby('bin', observed=False)[col].apply(lambda s: 100 * (s == 'conserved').mean())
        ax.plot(range(5), div.values, ls, color=HC_COLORS['divergent'], marker='o', ms=3, label=f'DS ({lab})')
        ax.plot(range(5), con.values, ls, color=HC_COLORS['conserved'], marker='o', ms=3, label=f'equivalent ({lab})')
    ax.set_xticks(range(5))
    ax.set_xticklabels(['<1', '1-5', '5-10', '10-20', '20-50'])
    ax.set_xlabel('minor isoform usage (%)')
    ax.set_ylabel('% of testable classifications\n(event x cell type)')
    ax.legend(fontsize=5, frameon=False, ncol=2, loc='lower center', bbox_to_anchor=(0.5, 1.0))
    label(ax, 'C', x=-0.25, y=1.2)

    pairs = [('H1_H2', 'H1 (GSE60424) PSI', 'H2 (GSE115736) PSI'), ('M1_M2', 'M1 PSI', 'M2 PSI'),
             ('H2_M1', 'Human (H2) PSI', 'Mouse (M1) PSI')]
    sgs = gs[2, :].subgridspec(1, 3, wspace=0.5)
    for i, (comp, xl, yl) in enumerate(pairs):
        ax = fig.add_subplot(sgs[i])
        j = jun[(jun['comparison'] == comp) & (jun['cell'] == 'NK') & jun['testable']]
        ax.hexbin(j['psi1'], j['psi2'], gridsize=30, bins='log', cmap='Greys', mincnt=1, linewidths=0)
        ax.plot([0, 1], [0, 1], color=KIND_COLORS[kind(comp)], lw=1)
        r = stats.pearsonr(j['psi1'], j['psi2'])[0]
        ax.set_title(f'{KIND_LABEL[kind(comp)]}, NK: r = {r:.2f}', fontsize=7)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_aspect('equal')
        if i == 0:
            label(ax, 'D', x=-0.35)
    fig.text(0.02, 0.0, 'Figure 3', fontsize=8)
    save(fig, 'Figure_3')


# ------------------------------------------------------------------ Figure 6
def figure6(ev, status):
    fib = ev[ev['comparison'] == 'FH_FM'].set_index('cluster')[STATUS]
    cls = status[status['status_high_low'] != 'not testable']
    fig = plt.figure(figsize=(7.2, 3.0))
    gs = fig.add_gridspec(1, 4, wspace=0.6)

    def bars(ax, rows, ns, colors, title, letter):
        df = pd.DataFrame(rows).reindex(columns=list(colors)).fillna(0)
        bottom = np.zeros(len(df))
        for s_, c in colors.items():
            ax.bar(df.index, df[s_], bottom=bottom, color=c, width=0.75, edgecolor='white', linewidth=0.5)
            bottom += df[s_].values
        for x, n in enumerate(ns):
            ax.text(x, 101, str(n), ha='center', va='bottom', fontsize=5)
        ax.set_ylim(0, 108)
        ax.set_title(title, fontsize=7)
        ax.set_ylabel('%')
        label(ax, letter, x=-0.35)

    def immune_to_fib(ax, klass, title, letter):  # status in fibroblasts of the events of an immune class
        rows, ns = [], []
        for cell in CELLS:
            f = cls.loc[(cls['cell'] == cell) & (cls['status_high_low'] == klass), 'cluster'].map(fib)
            f = f[f.notna() & (f != 'not_testable')]
            ns.append(len(f))
            rows.append(f.value_counts(normalize=True).rename(CELL_LABEL[cell]) * 100)
        bars(ax, rows, ns, STATUS_COLORS, title, letter)

    def fib_to_immune(ax, fstatus, title, letter):  # immune classification of the events of a fibroblast status
        clusters = fib[fib == fstatus].index
        rows, ns = [], []
        for cell in CELLS:
            c = cls[(cls['cell'] == cell) & cls['cluster'].isin(clusters)]['status_high_low']
            ns.append(len(c))
            rows.append(c.value_counts(normalize=True).rename(CELL_LABEL[cell]) * 100)
        bars(ax, rows, ns, HC_COLORS, title, letter)

    immune_to_fib(fig.add_subplot(gs[0, 0]), 'divergent', 'Immune DS:\nstatus in fibroblasts', 'A')
    immune_to_fib(fig.add_subplot(gs[0, 1]), 'conserved', 'Immune equivalent:\nstatus in fibroblasts', 'B')
    fib_to_immune(fig.add_subplot(gs[0, 2]), 'DS', 'Fibroblast DS:\nclass in immune cells', 'C')
    fib_to_immune(fig.add_subplot(gs[0, 3]), 'equivalent', 'Fibroblast equivalent:\nclass in immune cells', 'D')

    fig.legend([Patch(color=c) for c in STATUS_COLORS.values()], [f'Fibroblasts: {v}' for v in STATUS_LABEL.values()],
               loc='upper center', ncol=3, frameon=False, fontsize=6, bbox_to_anchor=(0.5, 1.12))
    fig.legend([Patch(color=c) for c in HC_COLORS.values()],
               ['Immune: differentially spliced', 'Immune: equivalent', 'Immune: inconclusive'],
               loc='upper center', ncol=3, frameon=False, fontsize=6, bbox_to_anchor=(0.5, 1.05))
    fig.text(0.02, -0.06, 'Figure 6', fontsize=8)
    save(fig, 'Figure_6')


def main():
    os.makedirs(OUT, exist_ok=True)
    np.random.seed(1)
    ev = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    ev['kind'] = ev['comparison'].map(kind)
    status = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
    jun = pd.read_csv(ED + 'effects/junction_effects.tsv', sep='\t')
    usage = pd.read_csv(ED + 'revision/minor_usage_by_event.tsv', sep='\t')
    n = figure2(ev, status, 'high', 'Figure_2', (3, 5))
    print('Figure 2G events:', n)
    n = figure2(ev, status, 'high_low', 'Figure_S2', (2, 2))
    print('Figure S2G events:', n)
    figure3(ev, jun, usage)
    figure6(ev, status)


if __name__ == '__main__':
    main()
