#!/usr/bin/env python3
"""
Within-species control comparisons for the reviewers' confounders (sex, age, tissue): each varies one factor
while cell type, lab and protocol are fixed. They are tested on the same AS events as the main analysis.

Steps (run in order):
  counts   add the junction counts of the control datasets (orthologous junction tables, human coordinates)
           to the events of the main analysis -> events/id50/AS_clusters_value_confounders.txt
           (the main analysis columns are unchanged; new columns are only added)
  design   groups files from the spec (confounder_comparisons.tsv) -> events/id50/leafcutter/confounders/
  summary  after LeafCutter and effect_size_analysis.py: compare each control comparison with the human-mouse and
           within-species comparisons of the same cell type on the same events -> events/id50/revision/confounders_*

Spec (ensembl116/metadata/confounder_comparisons.tsv): comparison, factor, species, cell, group, study, sample
(sample = column name in the counts table, i.e. '<sample>_<study>').
Usage: confounders_116.py counts <study>:<orthologous junction table> ...
       confounders_116.py design
       confounders_116.py summary
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
ED = B + 'events/id50/'
SPEC = B + 'metadata/confounder_comparisons.tsv'
COUNTS = ED + 'AS_clusters_value.txt'
COUNTS_EXT = ED + 'AS_clusters_value_confounders.txt'
LC = ED + 'leafcutter/confounders/'
CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']


def kind(c):
    if c.startswith('matched_'):
        return 'matched'
    if c == 'FH_FM':
        return 'fibro'
    if c.startswith(('human_', 'mouse_')):
        return 'control'
    a, b = c.split('_')
    return 'cross' if a[0] != b[0] else ('HH' if a[0] == 'H' else 'MM')


def counts(tables):
    from shared_junctions_116 import load_table
    base = pd.read_csv(COUNTS, sep='\t', index_col=0)
    key = base.index.str.rsplit(':', n=1).str[0]  # chr:start:end
    out = [base]
    for spec in tables:
        study, path = spec.split(':', 1)
        t = load_table(path)  # index chr:start:end, columns '<study>_<sample>'
        t.columns = [c.split('_', 1)[1] + '_' + c.split('_', 1)[0] for c in t.columns]  # -> '<sample>_<study>'
        add = t.reindex(key).fillna(0).astype(int)
        add.index = base.index
        print(study, f'{(add.sum(axis=1) > 0).sum()} of {len(add)} event junctions covered', flush=True)
        out.append(add)
    ext = pd.concat(out, axis=1)
    ext = ext.loc[:, ~ext.columns.duplicated()]
    ext.to_csv(COUNTS_EXT, sep='\t')
    print('wrote', COUNTS_EXT, ext.shape)


def design():
    spec = pd.read_csv(SPEC, sep='\t')
    cols = set(pd.read_csv(COUNTS_EXT, sep='\t', index_col=0, nrows=1).columns)
    missing = sorted(set(spec['sample']) - cols)
    if missing:
        sys.exit(f'samples not in {COUNTS_EXT}: {missing}')
    jobs = []
    for (comp, cell), g in spec.groupby(['comparison', 'cell']):
        if g['group'].nunique() != 2 or g.groupby('group').size().min() < 2:
            print('skip', comp, cell, g.groupby('group').size().to_dict())
            continue
        path = os.path.join(LC, comp, cell)
        os.makedirs(path, exist_ok=True)
        g[['sample', 'group']].to_csv(os.path.join(path, 'groups_file.txt'), sep='\t', header=False, index=False)
        jobs.append(path)
    open(os.path.join(LC, 'jobs.txt'), 'w').write('\n'.join(jobs) + '\n')
    print(len(jobs), 'LeafCutter runs; submit:  qsub -N lc_conf -t 1-%d run_leafcutter_116.sh %s %s'
          % (len(jobs), LC, COUNTS_EXT))


def summary():
    main = pd.read_csv(ED + 'effects/event_effects.tsv', sep='\t')
    ctrl = pd.read_csv(ED + 'effects_confounders/event_effects.tsv', sep='\t')
    used = set(pd.read_csv(SPEC, sep='\t')['comparison'])
    ctrl = ctrl[ctrl['comparison'].isin(used)]  # comparisons dropped from the spec (GSE124829) are not reported
    jm = pd.read_csv(ED + 'effects/junction_effects.tsv', sep='\t')
    jc = pd.read_csv(ED + 'effects_confounders/junction_effects.tsv', sep='\t')
    ev = pd.concat([main, ctrl], ignore_index=True)
    jun = pd.concat([jm, jc], ignore_index=True)
    ev['kind'] = ev['comparison'].map(kind)
    jun['kind'] = jun['comparison'].map(kind)
    rows = []
    for (comp, cell), c in ctrl.groupby(['comparison', 'cell']):
        tested = set(c.loc[c['testable'], 'cluster'])
        e = ev[(ev['cell'] == cell) & ev['kind'].isin(['control', 'cross', 'HH', 'MM']) & ev['cluster'].isin(tested)]
        e = e[e['comparison'].eq(comp) | ~e['kind'].eq('control')]
        # events testable in this control comparison and in every human-mouse pair of the cell type
        cross_ok = e[e['kind'].eq('cross') & e['testable']].groupby('cluster')['comparison'].nunique()
        common = cross_ok[cross_ok == cross_ok.max()].index if len(cross_ok) else []
        for k, g in e[e['cluster'].isin(common)].groupby('kind'):
            per = g.groupby('comparison').apply(lambda x: pd.Series({
                'pct_DS': 100 * (x['status_m0.05'] == 'DS').mean(),
                'median_max_abs_dpsi': x['max_abs_dpsi'].median()}))
            j = jun[(jun['cell'] == cell) & jun['cluster'].isin(common) & jun['testable'] & (jun['kind'] == k)]
            j = j[j['comparison'].eq(comp) | ~j['kind'].eq('control')]
            r = j.groupby('comparison').apply(lambda x: stats.pearsonr(x['psi1'], x['psi2'])[0]
                                              if len(x) > 2 else np.nan)
            rows.append(dict(control=comp, cell=cell, common_events=len(common), kind=k, comparisons=len(per),
                             pct_DS_median=per['pct_DS'].median(), pct_DS_min=per['pct_DS'].min(),
                             pct_DS_max=per['pct_DS'].max(),
                             median_max_abs_dpsi=per['median_max_abs_dpsi'].median(), pearson_median=r.median()))
        # paired test: events' |dPSI| in the control comparison vs the mean over human-mouse pairs
        p = e[e['cluster'].isin(common)].pivot_table(index='cluster', columns='comparison', values='max_abs_dpsi')
        cross_cols = [x for x in p if kind(x) == 'cross']
        if comp in p and cross_cols:
            q = p[[comp] + cross_cols].dropna()
            rows.append(dict(control=comp, cell=cell, common_events=len(q), kind='test',
                             wilcoxon_p_cross_gt_control=stats.wilcoxon(q[cross_cols].mean(1), q[comp],
                                                                        alternative='greater').pvalue))
    out = pd.DataFrame(rows)
    out.to_csv(ED + 'revision/confounders_summary.tsv', sep='\t', index=False)
    pd.set_option('display.width', 220)
    print(out.round(3).to_string(index=False))


if __name__ == '__main__':
    {'counts': lambda: counts(sys.argv[2:]), 'design': design, 'summary': summary}[sys.argv[1]]()
