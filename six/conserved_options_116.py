#!/usr/bin/env python3
"""
Options for the definition of a splicing conserved (unchanged) event, and the differentially spliced (DS) events of
the split definition, on the six-dataset events (events/id50_6ds).

Immune cells (per cell type; PSI of a junction in a sample = junction reads / event reads, samples with >= 20 reads;
testable = PSI in >= 2 samples of each of the six datasets; dataset PSI = mean over its samples):
  ANOVA: one-way ANOVA of the sample PSI of the junction with the largest range of dataset PSIs, six datasets as
         groups, Benjamini-Hochberg over the testable events of the cell type.
  DS (split definition): ANOVA adjusted P <= 0.05 and, in at least one junction of the event, gap >= 0.1 and
         1.5 x the largest within-species |dPSI| < the smallest human-mouse |dPSI|, where gap = |mean of the three
         human dataset PSIs - mean of the three mouse dataset PSIs| and the |dPSI| are between dataset PSIs.
         The same rule is applied to the nine splits of the datasets that mix the species (calibration).
  Conserved options, every junction of the event:
         range of the six dataset PSIs (max - min) <= 0.1 or <= 0.05;
         range of the PSI of all samples of the six datasets <= 0.1 or <= 0.05;
         each with ANOVA P (not adjusted) > 0.05 or > 0.1, or without the P value.
Single comparisons (fibroblasts, sex / age / tissue controls, and the within-species dataset pairs for reference;
LeafCutter p.adjust, max|dPSI| = largest |difference of the group mean PSI| over the junctions of the event):
  DS = p.adjust <= 0.05 and max|dPSI| >= 0.1;
  conserved options = max|dPSI| <= 0.05 or <= 0.1, with the LeafCutter P (not adjusted) > 0.05 or > 0.1, or without it;
  equivalent (90% CI within +-0.1) for reference.
Output: ensembl116/revision_6ds/conserved_options.xlsx (+ tsv)
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'separation'))
import build_supp_tables_116 as bst  # noqa: E402

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
E = B + 'events/id50_6ds/'
OUT = B + 'revision_6ds/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
HUMAN, MOUSE = ['H1', 'H2', 'H3'], ['M1', 'M2', 'M3']
DATASETS = HUMAN + MOUSE
MIN_READS, MIN_SAMPLES, FDR = 20, 2, 0.05
GAP, RATIO = 0.1, 1.5
SPLITS = [g for g in itertools.combinations(DATASETS, 3) if 'H1' in g]   # each split once; H1 H2 H3 = species


def cell_classes(counts, design, clu_all, cell):
    smp = design.index[design['cell'] == cell].intersection(counts.columns)
    groups = {d: design.index[(design['cell'] == cell) & (design['dataset'] == d)].intersection(smp) for d in DATASETS}
    tot = counts[smp].groupby(clu_all).transform('sum')
    psi = counts[smp] / tot.where(tot >= MIN_READS)
    ok = pd.concat([psi[g].notna().sum(axis=1) >= MIN_SAMPLES for g in groups.values()], axis=1).all(axis=1)
    ok = ok.groupby(clu_all).transform('all')
    psi, clu = psi[ok], clu_all[ok]
    allsmp = [s for g in groups.values() for s in g]
    means = pd.DataFrame({d: psi[g].mean(axis=1) for d, g in groups.items()})
    rng = means.max(axis=1) - means.min(axis=1)
    srng = psi[allsmp].max(axis=1) - psi[allsmp].min(axis=1)          # all samples of the six datasets
    junction = rng.groupby(clu).idxmax()
    p = {}
    for c, j in junction.items():
        v = [psi.loc[j, g].dropna().to_numpy() for g in groups.values()]
        allv = np.concatenate(v)
        p[c] = 1.0 if np.ptp(allv) == 0 else np.nan_to_num(f_oneway(*v).pvalue, nan=1.0)
    t = pd.DataFrame({'anova_p': pd.Series(p)})
    t['anova_padj'] = multipletests(t['anova_p'], method='fdr_bh')[1]
    t['dataset_range'] = rng.groupby(clu).max()
    t['sample_range'] = srng.groupby(clu).max()
    pairs = list(itertools.combinations(DATASETS, 2))
    for g1 in SPLITS:
        g2 = [d for d in DATASETS if d not in g1]
        gap = (means[list(g1)].mean(axis=1) - means[g2].mean(axis=1)).abs()
        win = pd.concat([(means[a] - means[b]).abs() for a, b in pairs if (a in g1) == (b in g1)], axis=1).max(axis=1)
        btw = pd.concat([(means[a] - means[b]).abs() for a, b in pairs if (a in g1) != (b in g1)], axis=1).min(axis=1)
        sep = ((gap >= GAP) & (RATIO * win < btw)).groupby(clu).any()
        t['split ' + ' '.join(g1) + ' vs ' + ' '.join(g2)] = sep & (t['anova_padj'] <= FDR)
    t['DS'] = t['split H1 H2 H3 vs M1 M2 M3']
    t['cell'] = cell
    return t.rename_axis('cluster').reset_index()


def main():
    counts = pd.read_csv(E + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(DATASETS)]
    clu_all = pd.Series(counts.index.str.split(':').str[3], index=counts.index)
    t = pd.concat([cell_classes(counts, design, clu_all, c) for c in CELLS], ignore_index=True)
    t.to_csv(OUT + 'conserved_options_events.tsv', sep='\t', index=False)
    opts = {}
    for what, col in (('dataset PSIs', 'dataset_range'), ('all sample PSIs', 'sample_range')):
        for thr in (0.1, 0.05):
            for pthr in (0.05, 0.1):
                opts[f'ANOVA P > {pthr} and range of {what} <= {thr}'] = (t['anova_p'] > pthr) & (t[col] <= thr)
            opts[f'range of {what} <= {thr} (no P value)'] = t[col] <= thr
    rows = []
    for cell, g in t.groupby('cell', sort=False):
        r = {'cell': cell, 'testable events': len(g), 'ANOVA P <= 0.05': int((g['anova_p'] <= 0.05).sum()),
             'ANOVA P > 0.05': int((g['anova_p'] > 0.05).sum()), 'ANOVA P > 0.1': int((g['anova_p'] > 0.1).sum()),
             'ANOVA adj P <= 0.05': int((g['anova_padj'] <= FDR).sum()),
             'DS (ANOVA adj P <= 0.05, gap >= 0.1, 1.5 x requirement)': int(g['DS'].sum())}
        for k, m in opts.items():
            n = int(m[g.index].sum())
            r[k] = f'{n} ({100 * n / len(g):.1f}%)'
        rows.append(r)
    allr = {'cell': 'events in >= 1 cell type', 'testable events': t['cluster'].nunique(),
            'ANOVA P <= 0.05': t.loc[t['anova_p'] <= 0.05, 'cluster'].nunique(),
            'ANOVA P > 0.05': t.loc[t['anova_p'] > 0.05, 'cluster'].nunique(),
            'ANOVA P > 0.1': t.loc[t['anova_p'] > 0.1, 'cluster'].nunique(),
            'ANOVA adj P <= 0.05': t.loc[t['anova_padj'] <= FDR, 'cluster'].nunique(),
            'DS (ANOVA adj P <= 0.05, gap >= 0.1, 1.5 x requirement)': t.loc[t['DS'], 'cluster'].nunique()}
    for k, m in opts.items():
        allr[k] = str(t.loc[m, 'cluster'].nunique())
    rows.append(allr)
    immune = pd.DataFrame(rows)
    # overlap of each conserved option with DS (should be none) and with ANOVA-significant events
    over = pd.DataFrame([{'option': k, 'classifications': int(m.sum()), 'also DS': int((m & t['DS']).sum()),
                          'ANOVA P <= 0.05': int((m & (t['anova_p'] <= 0.05)).sum())} for k, m in opts.items()])
    # DS calibration: the human-mouse split vs the nine splits that mix the species
    sp = [c for c in t.columns if c.startswith('split ')]
    cal = t.groupby('cell', sort=False)[sp].sum().T
    cal.index = cal.index.str.replace('split ', '')
    cal = cal.reset_index().rename(columns={'index': 'split'})
    mixed = cal[cal['split'] != 'H1 H2 H3 vs M1 M2 M3']
    cal.loc[len(cal)] = ['mean of the 9 mixed splits'] + [round(mixed[c].mean(), 1) for c in CELLS]
    # single comparisons
    ev = pd.concat([pd.read_csv(E + 'effects/event_effects.tsv', sep='\t'),
                    pd.read_csv(E + 'effects_confounders/event_effects.tsv', sep='\t')], ignore_index=True)
    used = set(pd.read_csv(B + 'metadata/confounder_comparisons.tsv', sep='\t')['comparison'])

    def kind(c):
        if c == 'FH_FM':
            return 'fibroblasts'
        if c in used:
            return 'control: ' + c.split('_')[1]
        a, b = c.split('_')[:2]
        return 'human-human' if a[0] == b[0] == 'H' else 'mouse-mouse' if a[0] == b[0] == 'M' else None
    ev['kind'] = ev['comparison'].map(lambda c: kind(c) if c.count('_') >= 1 else None)
    raw = []
    for sub in ('original', 'confounders'):
        for comp in os.listdir(E + f'leafcutter/{sub}'):
            for cell in os.listdir(E + f'leafcutter/{sub}/{comp}') if os.path.isdir(E + f'leafcutter/{sub}/{comp}') else []:
                f = E + f'leafcutter/{sub}/{comp}/{cell}/leafcutter_ds_cluster_significance.txt'
                if os.path.exists(f):
                    d = pd.read_csv(f, sep='\t', usecols=['cluster', 'p'])
                    raw.append(d.assign(cluster=d['cluster'].str.split(':').str[1], comparison=comp, cell=cell))
    raw = pd.concat(raw).rename(columns={'p': 'p_raw'})
    ev = ev.merge(raw, on=['cluster', 'comparison', 'cell'], how='left')
    ev = ev[ev['kind'].notna() & ev['testable']]
    sig = ev['p.adjust'] <= FDR
    # raw LeafCutter P; NaN (no LeafCutter P) is not counted in the P > options
    single = []
    for (k, comp, cell), g in ev.groupby(['kind', 'comparison', 'cell']):
        n = len(g)
        r = {'kind': k, 'comparison': comp, 'cell': cell, 'testable events': n,
             'DS (p.adjust <= 0.05, max|dPSI| >= 0.1)': int((sig[g.index] & (g['max_abs_dpsi'] >= 0.1)).sum()),
             'equivalent (90% CI within +-0.1)': int((g['status_m0.1'] == 'equivalent').sum())}
        for thr in (0.1, 0.05):
            for pthr in (0.05, 0.1):
                m1 = (g['p_raw'] > pthr) & (g['max_abs_dpsi'] <= thr)
                r[f'P > {pthr} and max|dPSI| <= {thr}'] = f'{int(m1.sum())} ({100 * m1.mean():.1f}%)'
            m2 = g['max_abs_dpsi'] <= thr
            r[f'max|dPSI| <= {thr} (no P value)'] = f'{int(m2.sum())} ({100 * m2.mean():.1f}%)'
        single.append(r)
    single = pd.DataFrame(single)
    order = {'fibroblasts': 0, 'control: sex': 1, 'control: age': 2, 'control: tissue': 3, 'human-human': 4,
             'mouse-mouse': 5}
    single = single.sort_values(['kind', 'cell'], key=lambda c: c.map(order) if c.name == 'kind' else c)
    sheets = {'immune_6datasets': immune, 'conserved_vs_DS_overlap': over, 'DS_split_calibration': cal,
              'single_comparisons': single}
    for k, v in sheets.items():
        v.to_csv(OUT + f'conserved_options_{k}.tsv', sep='\t', index=False)
    bst.OUT = OUT
    bst.write('conserved_options', 'Options for the definition of a splicing conserved event (six datasets)',
              'Immune cells: testable = PSI in >= 2 samples of each of the six datasets. ANOVA on the junction with the '
              'largest range of dataset PSIs (BH over the testable events). DS = ANOVA adj P <= 0.05 and, in at least '
              'one junction, gap (|mean of human dataset PSIs - mean of mouse dataset PSIs|) >= 0.1 and 1.5 x largest '
              'within-species |dPSI| < smallest human-mouse |dPSI| (dataset PSIs). Conserved options: range (max - min) '
              'of the dataset PSIs, or of the PSI of all samples, in every junction, <= 0.1 or <= 0.05, with ANOVA P (not '
              'adjusted) > 0.05 or > 0.1, or without the P value; n (% of testable). DS_split_calibration: the DS rule applied to every '
              'split of the six datasets into two groups of three. Single comparisons: LeafCutter p.adjust and '
              'max|dPSI| (group means); DS = p.adjust <= 0.05 and max|dPSI| >= 0.1; conserved options max|dPSI| <= '
              '0.1 or <= 0.05, with the LeafCutter P value (not adjusted) > 0.05 or > 0.1, or without the P value; events '
              'without a LeafCutter P value are not counted in the P options.', sheets)
    print(immune.to_string(index=False))
    print(cal.to_string(index=False))
    print(over.to_string(index=False))
    print(single[single['kind'] != 'human-human'].head(20).to_string(index=False))


if __name__ == '__main__':
    main()
