#!/usr/bin/env python3
"""
Analysis requested in DSA.docx (Tal Shay), to choose between the split analysis and the Wilcoxon analysis.

Per immune cell type, on the orthologous alternative splicing events (PSI of a junction in a sample = junction reads
/ event reads, samples with >= 20 reads in the event; an event is tested if it has a PSI in at least two samples of
every dataset):
  Junction: one junction per event, the one with the largest range of dataset PSIs (dataset PSI = mean PSI over the
            samples of the dataset; range = maximum minus minimum over the datasets; species labels not used).
  Stage 1 : one-way ANOVA of the junction PSI between the datasets; Benjamini-Hochberg adjusted P; events with
            adjusted P <= FDR (5% or 10%) go on.
  Stage 2 : all pairwise dataset comparisons: |dPSI| = |mean PSI of dataset a - mean PSI of dataset b| of the junction.
  Requirement: 2 x every within-species |dPSI| < every between-species |dPSI| (2 x max within < min between).
  Option 1: Wilcoxon rank-sum test of the junction PSI, all human vs all mouse samples, on the stage-1 events;
            Benjamini-Hochberg adjusted P <= FDR; with and without the requirement.
  Option 2: for each split of the datasets into a group with as many datasets as human and a group with as many as
            mouse (3-2: ten splits; the human-mouse split is one of them), for each event and each junction, the mean
            of the dataset PSIs of each group and the gap |group 1 - group 2|; an event is separated by the split if
            at least one junction has gap > T (T = 0.1, 0.2); with the requirement, that junction must also have
            2 x every within-group |dPSI| < every between-group |dPSI|. On the stage-1 events.
Output: ensembl116/revision_anova/dsa/DSA_tables_<n>datasets.xlsx (+ tsv per sheet)
Usage: dsa_analysis_116.py            (five datasets H1 H2 H3 M1 M2)
       dsa_analysis_116.py <counts> <design> <name> <dataset,...>  (e.g. with a third mouse dataset)
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'separation'))
import build_supp_tables_116 as bst  # noqa: E402

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
OUT = B + 'revision_anova/dsa/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
MIN_READS, MIN_SAMPLES = 20, 2
FDRS, GAPS = (0.05, 0.1), (0.1, 0.2)
OPT2_RATIO_ALT = 1.5   # option 2 with the requirement 1.5 x within-group < between-group
MIN_DPSI = 0.1   # option 1 requirement: the smallest between-species |dPSI| is also at least 0.1
EXAMPLE = ['clu_85337', 'PUF60']
EXAMPLE_JUNCTION = 'chr8:143818534:143821597'   # PUF60 exon 5 skipping junction (finds the event in other event sets)


def cell_tables(counts, design, clu_all, cell, datasets):
    human = [d for d in datasets if d.startswith('H')]
    mouse = [d for d in datasets if d.startswith('M')]
    smp = design.index[design['cell'] == cell].intersection(counts.columns)
    groups = {d: design.index[(design['cell'] == cell) & (design['dataset'] == d)].intersection(smp) for d in datasets}
    tot = counts[smp].groupby(clu_all).transform('sum')
    psi = counts[smp] / tot.where(tot >= MIN_READS)
    ok = pd.concat([psi[g].notna().sum(axis=1) >= MIN_SAMPLES for g in groups.values()], axis=1).all(axis=1)
    ok = ok.groupby(clu_all).transform('all')
    psi, clu = psi[ok], clu_all[ok]
    means = pd.DataFrame({d: psi[g].mean(axis=1) for d, g in groups.items()})
    drange = means.max(axis=1) - means.min(axis=1)
    junction = drange.groupby(clu).idxmax()                 # largest range of dataset means
    pairs = list(itertools.combinations(datasets, 2))
    splits = [g for g in itertools.combinations(datasets, len(human))
              if len(human) != len(mouse) or datasets[0] in g]          # 3 vs 3: each split once
    rows = []
    for c, j in junction.items():
        v = {d: psi.loc[j, g].dropna().to_numpy() for d, g in groups.items()}
        allv = np.concatenate(list(v.values()))
        p1 = 1.0 if np.ptp(allv) == 0 else np.nan_to_num(f_oneway(*v.values()).pvalue, nan=1.0)
        h = np.concatenate([v[d] for d in human]); m = np.concatenate([v[d] for d in mouse])
        pw = 1.0 if np.ptp(allv) == 0 else mannwhitneyu(h, m, alternative='two-sided').pvalue
        mu = means.loc[j]
        dp = {f'{a}-{b}': abs(mu[a] - mu[b]) for a, b in pairs}
        within = [dp[f'{a}-{b}'] for a, b in pairs if a[0] == b[0]]
        between = [dp[f'{a}-{b}'] for a, b in pairs if a[0] != b[0]]
        row = {'cluster': c, 'junction': j, 'anova_p': p1, 'wilcoxon_p': pw,
               'requirement': 2 * max(within) < min(between) and min(between) >= MIN_DPSI, **dp}
        for g1 in splits:   # option 1 for every split (human vs mouse is one of them): Wilcoxon + requirement
            a = np.concatenate([v[d] for d in g1]); b = np.concatenate([v[d] for d in datasets if d not in g1])
            row['p|' + ''.join(g1)] = 1.0 if np.ptp(allv) == 0 else mannwhitneyu(a, b, alternative='two-sided').pvalue
            wi = [dp[f'{x}-{y}'] for x, y in pairs if (x in g1) == (y in g1)]
            bt = [dp[f'{x}-{y}'] for x, y in pairs if (x in g1) != (y in g1)]
            row['req|' + ''.join(g1)] = 2 * max(wi) < min(bt) and min(bt) >= MIN_DPSI
        rows.append(row)
    t = pd.DataFrame(rows).set_index('cluster')
    t['anova_padj'] = multipletests(t['anova_p'], method='fdr_bh')[1]

    # option 2: all splits with the human / mouse group sizes; every junction of the event
    split_rows, split_rows15 = {}, {}
    for g1 in splits:
        g2 = [d for d in datasets if d not in g1]
        name = ' '.join(g1) + ' vs ' + ' '.join(g2)
        s1 = means[list(g1)].mean(axis=1)                    # mean of the dataset means of the group
        s2 = means[g2].mean(axis=1)
        gap = (s1 - s2).abs()
        win = [(a, b) for a, b in pairs if (a in g1) == (b in g1)]
        btw = [(a, b) for a, b in pairs if (a in g1) != (b in g1)]
        wmax = pd.concat([(means[a] - means[b]).abs() for a, b in win], axis=1).max(axis=1)
        bmin = pd.concat([(means[a] - means[b]).abs() for a, b in btw], axis=1).min(axis=1)
        req = 2 * wmax < bmin
        req15 = OPT2_RATIO_ALT * wmax < bmin                 # the same requirement with 1.5 instead of 2
        for T in GAPS:
            split_rows[(name, tuple(g1) == tuple(human), T, False)] = (gap > T).groupby(clu).any()
            split_rows15[(name, tuple(g1) == tuple(human), T, False)] = (gap > T).groupby(clu).any()
            split_rows15[(name, tuple(g1) == tuple(human), T, True)] = ((gap > T) & req15).groupby(clu).any()
            split_rows[(name, tuple(g1) == tuple(human), T, True)] = ((gap > T) & req).groupby(clu).any()
    return t, split_rows, psi, groups, means, junction, splits, split_rows15


def example_rows(cell, psi, groups, means, junction, datasets):
    c = EXAMPLE[0]
    js = psi.index[psi.index.str.endswith(':' + c)]
    rows = []
    for j in js:
        for d in datasets:
            v = psi.loc[j, groups[d]].dropna()
            rows.append({'cell': cell, 'junction': ':'.join(j.split(':')[:3]), 'chosen (largest dataset range)':
                         j == junction.get(c), 'dataset': d, 'sample PSIs': ', '.join(f'{x:.3f}' for x in v),
                         'dataset mean PSI': round(v.mean(), 3)})
    return rows


def main():
    if len(sys.argv) > 1:
        counts_path, design_path, name, datasets = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(',')
    else:
        counts_path, design_path = B + 'events/id50/AS_clusters_value.txt', B + 'events/id50/leafcutter/sample_design.tsv'
        name, datasets = '5datasets', ['H1', 'H2', 'H3', 'M1', 'M2']
    os.makedirs(OUT, exist_ok=True)
    counts = pd.read_csv(counts_path, sep='\t', index_col=0)
    design = pd.read_csv(design_path, sep='\t').set_index('sample')
    design = design[design['dataset'].isin(datasets)]
    clu_all = pd.Series(counts.index.str.split(':').str[3], index=counts.index)
    hit = [i for i in counts.index if i.startswith(EXAMPLE_JUNCTION + ':')]
    if hit:
        EXAMPLE[0] = hit[0].split(':')[3]
    s1, o1, o1s, o2, o215, ex, ex_calc = [], [], [], [], [], [], []
    for cell in CELLS:
        t, split_rows, psi, groups, means, junction, splits, split_rows15 = cell_tables(counts, design, clu_all, cell, datasets)
        s1.append({'cell': cell, 'tested events': len(t),
                   **{f'ANOVA adjusted P <= {f}': int((t['anova_padj'] <= f).sum()) for f in FDRS}})
        for f in FDRS:
            keep = t['anova_padj'] <= f
            w = t.loc[keep].copy()
            w['wilcoxon_padj'] = multipletests(w['wilcoxon_p'], method='fdr_bh')[1]
            sig = w['wilcoxon_padj'] <= f
            o1.append({'cell': cell, 'FDR': f, 'stage 1 events': int(keep.sum()), 'Wilcoxon significant': int(sig.sum()),
                       'Wilcoxon significant + requirement': int((sig & w['requirement']).sum())})
            for g1 in splits:
                sg = multipletests(w['p|' + ''.join(g1)].astype(float), method='fdr_bh')[1] <= f
                o1s.append({'cell': cell, 'FDR': f, 'split': ' '.join(g1) + ' vs ' + ' '.join(
                    d for d in datasets if d not in g1), 'human vs mouse': all(d.startswith('H') for d in g1),
                    'Wilcoxon significant': int(sg.sum()),
                    'Wilcoxon significant + requirement': int((sg & w['req|' + ''.join(g1)].astype(bool)).sum())})
            for (split, hm, T, req), sep in split_rows.items():
                if T == GAPS[0] and not req:
                    row = {'cell': cell, 'stage 1 FDR': f, 'split': split, 'human vs mouse': hm}
                    for T2 in GAPS:
                        for r2 in (False, True):
                            v = split_rows[(split, hm, T2, r2)].reindex(t.index[keep]).fillna(False)
                            row[f'gap > {T2}' + (' + requirement' if r2 else '')] = int(v.sum())
                    o2.append(row)
                    row = {'cell': cell, 'stage 1 FDR': f, 'split': split, 'human vs mouse': hm}
                    for T2 in GAPS:
                        for r2 in (False, True):
                            v = split_rows15[(split, hm, T2, r2)].reindex(t.index[keep]).fillna(False)
                            row[f'gap > {T2}' + (' + requirement' if r2 else '')] = int(v.sum())
                    o215.append(row)
            if f == 0.05 and EXAMPLE[0] in t.index:
                e = t.loc[EXAMPLE[0]]
                ww = w.loc[EXAMPLE[0]] if EXAMPLE[0] in w.index else None
                calc = {'cell': cell, 'chosen junction': ':'.join(e['junction'].split(':')[:3]),
                        'ANOVA P': e['anova_p'], 'ANOVA adjusted P': e['anova_padj'], 'Wilcoxon P': e['wilcoxon_p'],
                        'Wilcoxon adjusted P': None if ww is None else ww['wilcoxon_padj'],
                        'requirement (2 x max within < min between)': bool(e['requirement'])}
                calc.update({f'|dPSI| {k}': round(e[k], 3) for k in e.index if '-' in k and k[0] in 'HM'})
                for (split, hm, T, req), sep in split_rows.items():
                    calc[f'split {split}: gap > {T}' + (' + req' if req else '')] = bool(sep.get(EXAMPLE[0], False))
                calc['cell'] = cell
                ex_calc.append(calc)
                ex += example_rows(cell, psi, groups, means, junction, datasets)
        print(cell, 'done', flush=True)
    o1s = pd.DataFrame(o1s)   # one row per cell type and split: FDR 5% and 10%, with and without the requirement
    o1s = o1s.pivot_table(index=['cell', 'split', 'human vs mouse'], columns='FDR',
                          values=['Wilcoxon significant', 'Wilcoxon significant + requirement']).reset_index()
    o1s.columns = [a if not b else f'FDR {int(b * 100)}%: ' + a.replace('Wilcoxon significant', 'significant')
                   for a, b in o1s.columns]
    o1s = o1s[['cell', 'split', 'human vs mouse', 'FDR 5%: significant', 'FDR 5%: significant + requirement',
               'FDR 10%: significant', 'FDR 10%: significant + requirement']]
    o1s[o1s.columns[3:]] = o1s[o1s.columns[3:]].astype(int)
    o1s = o1s.sort_values(['cell', 'human vs mouse'], ascending=[True, False], key=lambda c: c.map(
        {k: i for i, k in enumerate(CELLS)}) if c.name == 'cell' else c)
    def option2(rows):
        o = pd.DataFrame(rows)
        o = o[o['stage 1 FDR'] == 0.05].drop(columns='stage 1 FDR')
        return o.sort_values(['cell', 'human vs mouse'], ascending=[True, False], key=lambda c: c.map(
            {k: i for i, k in enumerate(CELLS)}) if c.name == 'cell' else c)
    o2, o215 = option2(o2), option2(o215)
    exd = pd.DataFrame(ex)
    exd = exd[exd['chosen (largest dataset range)']].drop(columns='chosen (largest dataset range)')
    res = []
    for r in ex_calc:
        dps = {k: v for k, v in r.items() if k.startswith('|dPSI| ')}
        within = [v for k, v in dps.items() if k[7] == k[10]]
        between = [v for k, v in dps.items() if k[7] != k[10]]
        hm = ' '.join(d for d in datasets if d.startswith('H')) + ' vs ' + ' '.join(d for d in datasets
                                                                                   if d.startswith('M'))
        res.append({'cell': r['cell'], 'junction': r['chosen junction'],
                    'stage 1: ANOVA adjusted P': float(f"{r['ANOVA adjusted P']:.3g}"),
                    'option 1: Wilcoxon adjusted P': None if r['Wilcoxon adjusted P'] is None else
                    float(f"{r['Wilcoxon adjusted P']:.3g}"),
                    'max within-species |dPSI|': max(within), 'min between-species |dPSI|': min(between),
                    'requirement (2 x max within < min between)': r['requirement (2 x max within < min between)'],
                    'option 2 human vs mouse split: gap > 0.1': r[f'split {hm}: gap > 0.1'],
                    'gap > 0.2': r[f'split {hm}: gap > 0.2'],
                    'mixed splits with gap > 0.1 + requirement': sum(v for k, v in r.items() if k.startswith('split ')
                                                                    and hm not in k and k.endswith('gap > 0.1 + req'))})
    sheets = {'stage1_ANOVA': pd.DataFrame(s1), 'option1_Wilcoxon': pd.DataFrame(o1), 'option1_splits': o1s,
              'option2_splits': o2, 'option2_splits_1.5': o215,
              'example_PSI': exd, 'example_result': pd.DataFrame(res)}
    for k, v in sheets.items():
        v.to_csv(OUT + f'{k}_{name}.tsv', sep='\t', index=False)
    bst.OUT = OUT
    bst.write(f'DSA_tables_{name}', f'DSA analysis ({", ".join(datasets)})',
              'Junction: largest range of dataset mean PSIs. Stage 1: ANOVA between datasets, BH-adjusted P <= FDR. '
              'Option 1: Wilcoxon of the junction PSI, all human vs all mouse samples, BH-adjusted P <= FDR, on the '
              'stage-1 events; option1_splits: the same for every split of the datasets (Wilcoxon of the samples of '
              'one group vs the other, requirement within / between groups). Option 2: for each split of the datasets, mean of the dataset PSIs of each group, an '
              'event is counted if at least one junction has |gap| > T, on the stage-1 events. Requirement: 2 x every '
              'within-species (option 2: within-group) |dPSI| < every between-species (between-group) |dPSI|, dataset '
              f'means; in option 1 also the smallest between |dPSI| >= {MIN_DPSI}; option2_splits_1.5 / option_2_table_1.5: '
              f'option 2 with the requirement {OPT2_RATIO_ALT} x every within-group |dPSI| < every between-group |dPSI|. Option 2 at stage 1 FDR 5%. Example: {EXAMPLE[1]} exon 5 skipping, PSI of the chosen junction.', sheets)
    # wide tables (splits x cell types), one block per column of option2_splits; after the sheet of the long table
    import openpyxl
    from openpyxl.styles import Font
    x = OUT + f'DSA_tables_{name}.xlsx'
    wb = openpyxl.load_workbook(x)
    for long_name, wide_name, lab in (('option2_splits', 'option_2_table', 'requirement (2 x)'),
                                      ('option2_splits_1.5', 'option_2_table_1.5', 'requirement (1.5 x)')):
        o = sheets[long_name]
        ws = wb.create_sheet(wide_name, wb.sheetnames.index(long_name) + 1)
        blocks = [('gap > 0.1 + requirement', 'gap > 0.2 + requirement'), ('gap > 0.1', 'gap > 0.2')]
        r0 = 1
        for left, right in blocks:
            for c0, col in ((1, left), (8, right)):
                ws.cell(r0, c0, col.replace('requirement', lab)).font = Font(bold=True)
                ws.cell(r0 + 1, c0, 'split')
                for j, cell in enumerate(CELLS):
                    ws.cell(r0 + 1, c0 + 1 + j, cell)
                w = o.pivot(index='split', columns='cell', values=col)[CELLS]
                order = list(o.loc[o['human vs mouse'], 'split'].unique()) + \
                    sorted(o.loc[~o['human vs mouse'], 'split'].unique())
                for i, sp in enumerate(order):
                    ws.cell(r0 + 2 + i, c0, sp)
                    for j, cell in enumerate(CELLS):
                        ws.cell(r0 + 2 + i, c0 + 1 + j, int(w.loc[sp, cell]))
            r0 += len(order) + 3
        ws.column_dimensions['A'].width = ws.column_dimensions['H'].width = 22
    wb.save(x)
    bst.times_new_roman(x) if hasattr(bst, 'times_new_roman') else None
    for k in ('stage1_ANOVA', 'option1_Wilcoxon'):
        print(sheets[k].to_string(index=False))


if __name__ == '__main__':
    main()
