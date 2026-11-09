#!/usr/bin/env python3
"""
Additional analyses requested in the review, for one orthology variant.

  1 dpsi_cutoff_sensitivity.tsv    % DS per comparison type with max|dPSI| cutoffs 0.05 / 0.1 / 0.2 (R2 minor 3)
  2 common_events_summary.tsv      effect sizes on events testable in every comparison of a cell type (R2 1, 6)
  3 splice_site_subset_summary.tsv classification restricted to events with well-supported
                                   human-mouse splice-site correspondence (R4 3)
  4 coverage_bias_model.tsv        logistic regressions of testability / strict class on event coverage,
                                   number of junctions and events per gene (R2 minor 5)
  5 fibroblast_effects.tsv,        human-mouse dPSI in immune cell types vs fibroblasts on the same events, and a
    fibroblast_interaction.tsv     species x cell-type interaction test per event (R4 4)
  6 enrichment_<list>_<library>.tsv  MSigDB 2024.1.Hs enrichment of strict divergent-only / conserved-only genes
    enrichment_coverage_check.tsv     against all testable genes, and a coverage-adjusted check of significant terms
Usage: .venv/bin/python revision_analyses_116.py id50
Output: ensembl116/events/<variant>/revision/
"""
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']
MIN_COVERAGE = 20
MSIGDB_VERSION = '2024.1.Hs'
LIBRARIES = {'GO_Biological_Process_2025': 'c5.go.bp', 'Reactome_2022': 'c2.cp.reactome', 'hallmark': 'h.all'}
ENRICH_FDR = 0.1


def kind(c):
    if c.startswith('matched_'):
        return 'matched_' + kind(c[len('matched_'):])
    if c == 'FH_FM':
        return 'fibro'
    a, b = c.split('_')
    return 'cross' if a[0] != b[0] else ('HH' if a[0] == 'H' else 'MM')


def dpsi_cutoff_sensitivity(ev):
    rows = []
    t = ev[ev['testable'] | ev['lc_success']]
    for cut in (0.05, 0.1, 0.2):
        ds = t['lc_success'] & (t['p.adjust'] <= 0.05) & (t['max_abs_dpsi'] >= cut)
        g = t.assign(DS=ds).groupby(['kind', 'cell'])['DS'].mean().mul(100).rename(f'pct_DS_dpsi{cut}')
        rows.append(g)
    return pd.concat(rows, axis=1).round(1).reset_index()


def common_events_summary(ev, jun):
    rows = []
    e = ev[~ev['kind'].str.startswith('matched') & (ev['kind'] != 'fibro')]
    for cell in CELLS:
        ec = e[e['cell'] == cell]
        comps = ec['comparison'].unique()
        tested = ec[ec['testable']].groupby('cluster')['comparison'].nunique()
        common = tested[tested == len(comps)].index
        for comp in comps:
            x = ec[(ec['comparison'] == comp) & ec['cluster'].isin(common)]
            j = jun[(jun['comparison'] == comp) & (jun['cell'] == cell) & jun['cluster'].isin(common)]
            rows.append(dict(cell=cell, comparison=comp, kind=kind(comp), common_events=len(common),
                             pct_DS=100 * (x['status_m0.05'] == 'DS').mean(),
                             pct_equivalent_m005=100 * (x['status_m0.05'] == 'equivalent').mean(),
                             pct_equivalent_m01=100 * (x['status_m0.1'] == 'equivalent').mean(),
                             median_max_abs_dpsi=x['max_abs_dpsi'].median(),
                             pearson=stats.pearsonr(j['psi1'], j['psi2'])[0]))
    out = pd.DataFrame(rows)
    # paired test per cell: human-mouse event |dPSI| vs within-species, on the common events
    tests = []
    for cell in CELLS:
        p = e[(e['cell'] == cell)].pivot_table(index='cluster', columns='comparison', values='max_abs_dpsi')
        p = p.dropna()
        cross = [c for c in p if kind(c) == 'cross']
        within = [c for c in p if kind(c) in ('HH', 'MM')]
        tests.append(dict(cell=cell, events=len(p),
                          wilcoxon_p_cross_gt_within=stats.wilcoxon(p[cross].mean(1), p[within].mean(1),
                                                                    alternative='greater').pvalue))
    return out.round(3), pd.DataFrame(tests)


def splice_site_subset(status):
    rows = []
    for subset, s in (('all', status), ('correspondence', status[status['splice_site_correspondence'] == True]),
                      ('no_correspondence', status[status['splice_site_correspondence'] == False])):
        for cell, g in s.groupby('cell'):
            t = g[g['status_high'] != 'not testable']
            rows.append(dict(subset=subset, cell=cell, testable=len(t),
                             pct_high_divergent=100 * (t['status_high'] == 'divergent').mean(),
                             pct_high_conserved=100 * (t['status_high'] == 'conserved').mean(),
                             pct_high_low_divergent=100 * (t['status_high_low'] == 'divergent').mean(),
                             pct_high_low_conserved=100 * (t['status_high_low'] == 'conserved').mean(),
                             pct_cross_DS=100 * (t['cross_DS'] / t['cross_testable']).mean(),
                             pct_within_DS=100 * (t['within_DS'] / t['within_testable'].replace(0, np.nan)).mean()))
    return pd.DataFrame(rows).round(2)


def event_coverage(counts):
    """Mean reads per sample of each event (sum of its junctions), over all immune samples."""
    immune = [c for c in counts.columns if 'EMTAB5919' not in c]
    per_cluster = counts[immune].groupby(counts.index.str.rsplit(':', n=1).str[1]).sum()
    return np.log10(per_cluster.mean(axis=1) + 1).rename('log10_coverage')


def coverage_bias_model(status, ann, cov):
    d = status.merge(ann[['cluster', 'n_junctions']], on='cluster', how='left')
    d = d.merge(cov, left_on='cluster', right_index=True, how='left')
    d['events_per_gene'] = d['ensembl'].map(ann.groupby('ensembl').size())
    d['testable'] = (d['status_high'] != 'not testable').astype(int)
    d['conserved'] = (d['status_high'] == 'conserved').astype(int)
    d['divergent'] = (d['status_high'] == 'divergent').astype(int)
    rows = []
    for outcome, data in (('testable', d), ('conserved', d[d['testable'] == 1]), ('divergent', d[d['testable'] == 1])):
        m = smf.logit(f'{outcome} ~ log10_coverage + n_junctions + events_per_gene + C(cell)', data=data).fit(disp=0)
        for term in ('log10_coverage', 'n_junctions', 'events_per_gene'):
            lo, hi = np.exp(m.conf_int().loc[term])
            rows.append(dict(outcome=outcome, n=int(m.nobs), term=term, odds_ratio=np.exp(m.params[term]),
                             ci_low=lo, ci_high=hi, p=m.pvalues[term]))
    return pd.DataFrame(rows), d


def group_psi(counts, samples):
    sub = counts[samples]
    clusters = counts.index.str.rsplit(':', n=1).str[1]
    totals = sub.groupby(clusters).transform('sum')
    return sub / totals.where(totals >= MIN_COVERAGE)


def fibroblast_analyses(counts, lc_dir, ev):
    """Human-mouse dPSI in fibroblasts vs each immune cell type (H2-M1, as in the manuscript), same events."""
    groups = lambda path: pd.read_csv(path, sep='\t', header=None, names=['s', 'g'])
    gf = groups(os.path.join(lc_dir, 'original', 'FH_FM', 'Fibro', 'groups_file.txt'))
    fh, fm = gf.loc[gf['g'] == 'FH', 's'].tolist(), gf.loc[gf['g'] == 'FM', 's'].tolist()
    clusters = counts.index.str.rsplit(':', n=1).str[1]
    eff, inter = [], []
    for cell in CELLS:
        gi = groups(os.path.join(lc_dir, 'original', 'H2_M1', cell, 'groups_file.txt'))
        ih, im = gi.loc[gi['g'] == 'H2', 's'].tolist(), gi.loc[gi['g'] == 'M1', 's'].tolist()
        psi = group_psi(counts, ih + im + fh + fm)
        stats_ = {}
        for name, s in (('ih', ih), ('im', im), ('fh', fh), ('fm', fm)):
            x = psi[s]
            stats_[name] = (x.mean(1), x.var(1, ddof=1), x.notna().sum(1))
        d_imm = stats_['ih'][0] - stats_['im'][0]
        d_fib = stats_['fh'][0] - stats_['fm'][0]
        v = sum(stats_[k][1] / stats_[k][2] for k in stats_)
        df = v ** 2 / sum((stats_[k][1] / stats_[k][2]) ** 2 / (stats_[k][2] - 1) for k in stats_)
        ok = pd.concat([stats_[k][2] >= 2 for k in stats_], axis=1).all(1) & (v > 0)
        t = (d_imm - d_fib) / np.sqrt(v)
        p = 2 * stats.t.sf(np.abs(t), df)
        j = pd.DataFrame({'cluster': clusters, 'd_immune': d_imm, 'd_fibro': d_fib, 'p': p}).loc[ok]
        # event level: junction with the largest |interaction|, Bonferroni within the event
        j['abs_int'] = (j['d_immune'] - j['d_fibro']).abs()
        n_j = j.groupby('cluster')['cluster'].transform('size')
        j['p_event'] = np.minimum(j['p'] * n_j, 1)
        e = j.sort_values('abs_int', ascending=False).groupby('cluster').first().reset_index()
        e['fdr'] = multipletests(e['p_event'], method='fdr_bh')[1]
        e['cell'] = cell
        e['interaction'] = (e['fdr'] <= 0.05) & (e['abs_int'] >= 0.1)
        st = ev[(ev['comparison'] == 'H2_M1') & (ev['cell'] == cell)].set_index('cluster')['status_m0.1']
        fst = ev[(ev['comparison'] == 'FH_FM')].set_index('cluster')['status_m0.1']
        e['status_immune_H2_M1'] = e['cluster'].map(st)
        e['status_fibro'] = e['cluster'].map(fst)
        inter.append(e)
        eff.append(dict(cell=cell, events=len(e),
                        pearson_dpsi_immune_vs_fibro=stats.pearsonr(j['d_immune'], j['d_fibro'])[0],
                        median_event_abs_dpsi_immune=e['d_immune'].abs().median(),
                        median_event_abs_dpsi_fibro=e['d_fibro'].abs().median(),
                        wilcoxon_p_immune_gt_fibro=stats.wilcoxon(e['d_immune'].abs(), e['d_fibro'].abs(),
                                                                  alternative='greater').pvalue,
                        pct_interaction=100 * e['interaction'].mean(),
                        pct_interaction_among_immune_DS=100 * e.loc[e['status_immune_H2_M1'] == 'DS',
                                                                     'interaction'].mean()))
    return pd.DataFrame(eff), pd.concat(inter, ignore_index=True)


def enrichment(status, ann, cov, out_dir):
    import gseapy as gp
    from gseapy import Msigdb
    s = status[status['symbol'].fillna('') != '']
    testable = s[s['status_high'] != 'not testable']
    background = sorted(set(testable['symbol']))
    div = set(testable.loc[testable['status_high'] == 'divergent', 'symbol'])
    con = set(testable.loc[testable['status_high'] == 'conserved', 'symbol'])
    lists = {'divergent_only': sorted(div - con), 'conserved_only': sorted(con - div)}
    msig = Msigdb()
    summary, checks = [], []
    gene_cov = s.merge(cov, left_on='cluster', right_index=True).groupby('symbol')['log10_coverage'].max()
    for lib, cat in LIBRARIES.items():
        gmt = msig.get_gmt(category=cat, dbver=MSIGDB_VERSION)
        for name, genes in lists.items():
            res = gp.enrich(gene_list=genes, gene_sets=gmt, background=background, outdir=None).results
            res = res.sort_values('Adjusted P-value')
            res.to_csv(os.path.join(out_dir, f'enrichment_{name}_{lib}.tsv'), sep='\t', index=False)
            sig = res[res['Adjusted P-value'] <= ENRICH_FDR]
            summary.append(dict(list=name, genes=len(genes), background=len(background), library=lib,
                                significant_terms=len(sig)))
            # coverage-adjusted check: is term membership still associated with the class
            # after adjusting for the gene's event coverage?
            for _, r in sig.head(20).iterrows():
                members = set(gmt[r['Term']])
                d = pd.DataFrame({'symbol': background})
                d['in_list'] = d['symbol'].isin(genes).astype(int)
                d['member'] = d['symbol'].isin(members).astype(int)
                d['log10_coverage'] = d['symbol'].map(gene_cov)
                m = smf.logit('in_list ~ member + log10_coverage', data=d.dropna()).fit(disp=0)
                checks.append(dict(list=name, library=lib, term=r['Term'], fdr=r['Adjusted P-value'],
                                   overlap=r['Overlap'], adjusted_odds_ratio=np.exp(m.params['member']),
                                   adjusted_p=m.pvalues['member']))
    pd.DataFrame(checks).to_csv(os.path.join(out_dir, 'enrichment_coverage_check.tsv'), sep='\t', index=False)
    genes_mixed = len(div & con)
    return pd.DataFrame(summary), genes_mixed, lists


def main():
    variant = sys.argv[1]
    ed = BASE + f'events/{variant}/'
    out = ed + 'revision/'
    os.makedirs(out, exist_ok=True)
    ev = pd.read_csv(ed + 'effects/event_effects.tsv', sep='\t')
    ev['kind'] = ev['comparison'].map(kind)
    jun = pd.read_csv(ed + 'effects/junction_effects.tsv', sep='\t')
    status = pd.read_csv(ed + 'classification/event_status_by_cell.tsv', sep='\t')
    ann = pd.read_csv(ed + 'classification/events_annotation.tsv', sep='\t')
    counts = pd.read_csv(ed + 'AS_clusters_value.txt', sep='\t', index_col=0)
    cov = event_coverage(counts)

    r1 = dpsi_cutoff_sensitivity(ev)
    r1.to_csv(out + 'dpsi_cutoff_sensitivity.tsv', sep='\t', index=False)
    print(r1[~r1['kind'].str.startswith('matched')].to_string(index=False))

    r2, t2 = common_events_summary(ev, jun)
    r2.to_csv(out + 'common_events_summary.tsv', sep='\t', index=False)
    t2.to_csv(out + 'common_events_tests.tsv', sep='\t', index=False)
    print(r2.groupby(['kind', 'cell'])[['common_events', 'pct_DS', 'pct_equivalent_m005', 'median_max_abs_dpsi',
                                        'pearson']].agg(['min', 'max']).round(2).to_string())
    print(t2.to_string(index=False))

    r3 = splice_site_subset(status)
    r3.to_csv(out + 'splice_site_subset_summary.tsv', sep='\t', index=False)
    print(r3.to_string(index=False))

    r4, d4 = coverage_bias_model(status, ann, cov)
    r4.to_csv(out + 'coverage_bias_model.tsv', sep='\t', index=False)
    print(r4.round(4).to_string(index=False))

    r5, i5 = fibroblast_analyses(counts, ed + 'leafcutter', ev)
    r5.to_csv(out + 'fibroblast_effects.tsv', sep='\t', index=False)
    i5.to_csv(out + 'fibroblast_interaction.tsv', sep='\t', index=False)
    print(r5.round(4).to_string(index=False))

    r6, mixed, lists = enrichment(status, ann, cov, out)
    r6.to_csv(out + 'enrichment_summary.tsv', sep='\t', index=False)
    print(r6.to_string(index=False))
    print('genes with both strict divergent and strict conserved events:', mixed)


if __name__ == '__main__':
    main()
