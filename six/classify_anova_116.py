#!/usr/bin/env python3
"""
Human-mouse classification of each orthologous alternative splicing event in each immune cell type, three human
(H1-H3) and three mouse (M1-M3) datasets (DSA option 1 with the requirement).

ANOVA junction: for each event and cell type, the junction with the largest range of dataset PSIs (maximum minus
minimum over the six datasets; the species labels are not used). PSI of a junction in a sample = junction reads /
reads of the event (LeafCutter cluster), in samples with >= 20 reads in the event; dataset PSI = mean over those samples.
Testable: the event has a PSI in at least two samples of each of the six datasets.
ANOVA: LeafCutter Dirichlet-multinomial model of all junctions of the event, six datasets as six groups, likelihood-
  ratio test (leafcutter_anova6_116.R -> events/<variant>_anova/leafcutter_anova6/<cell>.tsv); P adjust =
  Benjamini-Hochberg over the testable events of the cell type; events that the LeafCutter ANOVA did not test
  (LeafCutter filters) are not testable. 'junction' = the junction with the largest range of dataset PSIs (reported, not tested).
Splicing conserved ('conserved'): ANOVA P (not adjusted) > 0.1 and, in every junction, the dataset PSIs within a range
  of at most 0.05.
Differentially spliced ('divergent'): ANOVA adjusted P <= 0.05 and, in at least one junction, gap >= 0.1 and
  1.5 x the largest within-species |dPSI| < the smallest human-mouse |dPSI| (gap = |mean of the human dataset PSIs -
  mean of the mouse dataset PSIs|; |dPSI| between dataset PSIs: 3 human-human, 3 mouse-mouse, 9 human-mouse pairs).
All other testable events: 'inconclusive'.
status_high, status_high_low and status_confidence carry the same class (downstream scripts are shared).
Output: events/<variant>_anova/classification/{events_annotation, event_status_by_cell, summary_by_cell,
        event_tests}.tsv
Usage: classify_anova_116.py <variant> (id50_6ds, id70_6ds, id80_6ds, id50_rbh_6ds)
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import classify_events_116 as ce  # noqa: E402

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
HUMAN, MOUSE = ['H1', 'H2', 'H3'], ['M1', 'M2', 'M3']
DATASETS = HUMAN + MOUSE
MIN_READS, MIN_SAMPLES = 20, 2
FDR, RANGE, RATIO, MIN_GAP, CONS_P = 0.05, 0.05, 1.5, 0.1, 0.1
PAIRS = list(itertools.combinations(DATASETS, 2))
HM_PAIRS = [p for p in PAIRS if p[0][0] != p[1][0]]
WITHIN_PAIRS = [p for p in PAIRS if p[0][0] == p[1][0]]
LINKS = ['AS_clusters_value.txt', 'AS_clusters_value_confounders.txt', 'junction_presence.tsv', 'effects',
         'effects_confounders', 'leafcutter', 'sample_design.tsv']


LC_P = {}   # cell -> {cluster: LeafCutter P} (filled in main)


def lc_pvalues(variant, cell):
    t = pd.read_csv(BASE + f'events/{variant}_anova/leafcutter_anova6/{cell}.tsv', sep='\t')
    return dict(zip(t['cluster'], t['p']))


def classify_cell(counts, design, clu_all, cell):
    smp = design.index[design['cell'] == cell].intersection(counts.columns)
    groups = {d: design.index[(design['cell'] == cell) & (design['dataset'] == d)].intersection(smp)
              for d in DATASETS}
    tot = counts[smp].groupby(clu_all).transform('sum')
    psi = counts[smp] / tot.where(tot >= MIN_READS)
    n_ok = pd.DataFrame({d: psi[g].notna().sum(axis=1) for d, g in groups.items()})
    testable_j = (n_ok >= MIN_SAMPLES).all(axis=1).groupby(clu_all).transform('all')
    psi = psi[testable_j]
    clu = clu_all[testable_j]
    means = pd.DataFrame({d: psi[g].mean(axis=1) for d, g in groups.items()})
    rng = means.max(axis=1) - means.min(axis=1)
    junction = rng.groupby(clu).idxmax()                           # ANOVA junction: largest range
    max_range = rng.groupby(clu).max()
    # split rule per junction: gap between the mean of the human and of the mouse dataset PSIs, and the dataset |dPSI|
    gap = (means[HUMAN].mean(axis=1) - means[MOUSE].mean(axis=1)).abs()
    wmax = pd.concat([(means[a] - means[b]).abs() for a, b in WITHIN_PAIRS], axis=1).max(axis=1)
    bmin = pd.concat([(means[a] - means[b]).abs() for a, b in HM_PAIRS], axis=1).min(axis=1)
    jsep = (gap >= MIN_GAP) & (RATIO * wmax < bmin)
    # junction reported for the species difference: a separating junction with the largest gap, else the largest gap
    score = gap + jsep.astype(float) * 10
    ds_junction = score.groupby(clu).idxmax()
    rows = []
    for c, j in junction.items():
        p_anova = LC_P[cell].get(c, np.nan)          # LeafCutter, six groups (NaN: not tested by LeafCutter)
        mu = means.loc[j]
        k = ds_junction[c]
        rows.append({'cell': cell, 'cluster': c, 'junction': j, 'n_junctions': int((clu == c).sum()),
                     **{f'PSI_{d}': round(mu[d], 4) for d in DATASETS},
                     'max_range': max_range[c], 'anova_p': p_anova,
                     'species_junction': k, 'gap': gap[k], 'max_within_dpsi': wmax[k], 'min_HM_dpsi': bmin[k],
                     'separated': bool(jsep[k])})
    t = pd.DataFrame(rows)
    ok = t['anova_p'].notna()
    t['anova_p_adjust'] = np.nan
    t.loc[ok, 'anova_p_adjust'] = multipletests(t.loc[ok, 'anova_p'], method='fdr_bh')[1]
    t['unchanged'] = (t['anova_p'] > CONS_P) & (t['max_range'] <= RANGE)
    t['DS'] = (t['anova_p_adjust'] <= FDR) & t['separated']
    t['status'] = np.select([t['DS'], t['unchanged']], ['divergent', 'conserved'], 'inconclusive')
    t.loc[t['anova_p'].isna(), 'status'] = 'not testable'   # not tested by the LeafCutter ANOVA: untested, not inconclusive
    return t


def main():
    variant = sys.argv[1]
    new = BASE + f'events/{variant}_anova/'
    os.makedirs(new + 'classification', exist_ok=True)
    os.makedirs(new + 'revision', exist_ok=True)
    for name in LINKS:
        if os.path.exists(BASE + f'events/{variant}/{name}') and not os.path.lexists(new + name):
            os.symlink(BASE + f'events/{variant}/{name}', new + name)
    old = BASE + f'events/{variant}/classification/'
    if os.path.exists(old + 'events_annotation.tsv'):
        ann = pd.read_csv(old + 'events_annotation.tsv', sep='\t')
    else:
        ann = ce.annotate(variant)
    ann.to_csv(new + 'classification/events_annotation.tsv', sep='\t', index=False)
    counts = pd.read_csv(BASE + f'events/{variant}/AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(BASE + f'events/{variant}/sample_design.tsv', sep='\t').set_index('sample')
    design = design[design['dataset'].isin(DATASETS)]
    clu_all = pd.Series(counts.index.str.split(':').str[3], index=counts.index)
    for c in CELLS:
        LC_P[c] = lc_pvalues(variant, c)
    tests = pd.concat([classify_cell(counts, design, clu_all, c) for c in CELLS], ignore_index=True)
    tests.to_csv(new + 'classification/event_tests.tsv', sep='\t', index=False, float_format='%.6g')
    if os.path.exists(old + 'event_status_by_cell.tsv'):   # per-pair statuses (LeafCutter), id50 only
        status = pd.read_csv(old + 'event_status_by_cell.tsv', sep='\t')
        status = status.drop(columns=[c for c in ('status_confidence', 'status_high', 'status_high_low')
                                      if c in status.columns])
    else:
        status = pd.MultiIndex.from_product([CELLS, ann['cluster']], names=['cell', 'cluster']).to_frame(index=False)
        status = status.merge(ann[['cluster', 'symbol', 'ensembl', 'event_type', 'splice_site_correspondence']],
                              on='cluster', how='left')
        status['DS_within_species'] = False
    status = status.merge(tests[['cell', 'cluster', 'status', 'junction', 'anova_p', 'anova_p_adjust',
                                 'species_junction', 'gap', 'max_within_dpsi', 'min_HM_dpsi', 'max_range']],
                          on=['cell', 'cluster'], how='left')
    status['status'] = status['status'].fillna('not testable')
    for col in ('status_confidence', 'status_high', 'status_high_low'):
        status[col] = status['status']
    status = status.drop(columns='status')
    status.to_csv(new + 'classification/event_status_by_cell.tsv', sep='\t', index=False)
    summary = pd.concat([status.groupby('cell')[f'status_{n}'].value_counts().unstack(fill_value=0)
                         .add_prefix(f'{n}_') for n in ('high', 'high_low')], axis=1)
    for n in ('divergent', 'conserved'):
        summary[f'high_low_{n}_DS_within_species'] = (status['DS_within_species'].fillna(False).astype(bool) &
                                                     (status['status_high_low'] == n)).groupby(status['cell']).sum()
    summary.to_csv(new + 'classification/summary_by_cell.tsv', sep='\t')
    print(variant)
    print(summary[[c for c in summary.columns if c.startswith('high_') and not c.startswith('high_low')]].to_string())
    print('unique events: DS', status.loc[status['status_high'] == 'divergent', 'cluster'].nunique(),
          'unchanged', status.loc[status['status_high'] == 'conserved', 'cluster'].nunique(),
          'testable', status.loc[status['status_high'] != 'not testable', 'cluster'].nunique())



if __name__ == '__main__':
    main()
