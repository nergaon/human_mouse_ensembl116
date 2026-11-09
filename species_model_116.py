#!/usr/bin/env python3
"""
Species effect tested against the variation between datasets of the same species (reviewer 2: separate the
species effect from within-species dataset variability).

Model, per junction and cell type: PSI(sample) = mu + beta * species + u(dataset) + e, with dataset as a random
effect. Species is constant within a dataset, so beta is estimated and tested at the dataset level (the two-stage
/ summary-statistics form of this model): the mean PSI of each dataset (samples with >= 20 reads in the event)
is one observation, and
    beta  = mean of the human dataset means - mean of the mouse dataset means
    SE    = s * sqrt(1/nH + 1/nM),  s^2 = pooled variance of the dataset means within each species
    df    = nH + nM - 2
so the uncertainty of beta includes both the dataset effects and the replicate noise.
A dataset is used when it has at least 2 such samples (and 60% of its samples, the LeafCutter -g rule); an event
is analysed in a cell type when >= 2 human datasets and both mouse datasets are usable for all its junctions.

Event status (per cell type):
  divergent   species effect significant (event p = min junction p x number of junctions, BH across the events
              of the cell type, FDR <= 0.05) and max|beta| >= 0.1
  conserved   the 90% CI of beta of every junction within +-0.1 (TOST at 0.05); also reported for +-0.05
  inconclusive  otherwise
Usage: species_model_116.py <variant>
Output: events/<variant>/revision/species_model_junctions.tsv, species_model_events.tsv, species_model_summary.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import effect_size_analysis as esa  # noqa: E402

BASE = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/'
HUMAN, MOUSE = ['H1', 'H2', 'H3'], ['M1', 'M2']
CELL = {'BCell': 'NveB', 'Cd8T': 'CD8T'}
FDR, DPSI, MARGINS = 0.05, 0.1, (0.05, 0.1)


def dataset_samples(ld):
    """{(cell, dataset): samples} from the original LeafCutter groups files."""
    out = {}
    for comp in os.listdir(ld):
        if comp == 'FH_FM':
            continue
        for c in os.listdir(os.path.join(ld, comp)):
            gf = os.path.join(ld, comp, c, 'groups_file.txt')
            if os.path.exists(gf):
                g = pd.read_csv(gf, sep='\t', header=None)
                for lab, x in g.groupby(1):
                    out[(CELL.get(c, c), lab)] = x[0].tolist()
    return out


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def analyse_cell(counts, samples, cell):
    means = {}
    for ds in HUMAN + MOUSE:
        s = samples[(cell, ds)]
        sub = counts[s + ['cluster']]
        tot = sub.groupby('cluster')[s].transform('sum')
        psi = sub[s] / tot.where(tot >= esa.MIN_COVERAGE)
        n = psi.notna().sum(1)
        need = max(2, int(np.floor(0.6 * len(s))))
        means[ds] = psi.mean(1).where(n >= need)
    m = pd.DataFrame(means)
    m['cluster'] = counts['cluster']
    # a dataset is used for an event only if usable for all junctions of the event
    usable = m[HUMAN + MOUSE].notna().groupby(m['cluster']).transform('all')
    m[HUMAN + MOUSE] = m[HUMAN + MOUSE].where(usable)
    h, mo = m[HUMAN], m[MOUSE]
    nh, nm = h.notna().sum(1), mo.notna().sum(1)
    beta = h.mean(1) - mo.mean(1)
    ss = h.var(1, ddof=1).fillna(0) * (nh - 1) + mo.var(1, ddof=1).fillna(0) * (nm - 1)
    df = nh + nm - 2
    s2 = ss / df
    se = np.sqrt(s2 * (1 / nh + 1 / nm))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), df)
    p = pd.Series(p, index=m.index).where(se > 0, np.where(beta.abs() > 0, 0.0, 1.0))
    half = stats.t.ppf(0.95, df) * se
    jun = pd.DataFrame({'cell': cell, 'cluster': m['cluster'], 'n_human_datasets': nh, 'n_mouse_datasets': nm,
                        'beta': beta, 'se': se, 'df': df, 'p': p, 'ci_low': beta - half, 'ci_high': beta + half,
                        'sd_within_species': np.sqrt(s2)})
    for ds in HUMAN + MOUSE:
        jun[f'mean_{ds}'] = m[ds]
    jun = jun[(nh >= 2) & (nm == 2)]
    jun['ci_extreme'] = np.maximum(jun['ci_low'].abs(), jun['ci_high'].abs())
    ev = jun.groupby('cluster').agg(cell=('cell', 'first'), n_junctions=('beta', 'size'),
                                    n_human_datasets=('n_human_datasets', 'min'),
                                    max_abs_beta=('beta', lambda x: x.abs().max()), min_p=('p', 'min'),
                                    max_ci=('ci_extreme', 'max'), max_sd_within=('sd_within_species', 'max'))
    ev['p_event'] = np.minimum(ev['min_p'] * ev['n_junctions'], 1)
    ev['fdr'] = bh(ev['p_event'])
    div = (ev['fdr'] <= FDR) & (ev['max_abs_beta'] >= DPSI)
    for mg in MARGINS:
        ev[f'status_m{mg}'] = np.select([div, ev['max_ci'] < mg], ['divergent', 'conserved'], 'inconclusive')
    return jun.reset_index(), ev.reset_index()


def main():
    variant = sys.argv[1]
    E = BASE + variant + '/'
    counts = esa.read_counts(E + 'AS_clusters_value.txt')
    samples = dataset_samples(E + 'leafcutter/original')
    juns, evs = [], []
    for cell in ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']:
        j, e = analyse_cell(counts, samples, cell)
        juns.append(j)
        evs.append(e)
        print(cell, len(e), e['status_m0.1'].value_counts().to_dict(), flush=True)
    jun, ev = pd.concat(juns), pd.concat(evs)
    out = E + 'revision/'
    jun.to_csv(out + 'species_model_junctions.tsv', sep='\t', index=False)
    ev.to_csv(out + 'species_model_events.tsv', sep='\t', index=False)
    summ = pd.concat([ev.groupby('cell')[f'status_m{m}'].value_counts().unstack(fill_value=0).add_suffix(f'_m{m}')
                      for m in MARGINS], axis=1)
    summ['events'] = ev.groupby('cell').size()
    summ.to_csv(out + 'species_model_summary.tsv', sep='\t')
    print(summ.to_string())


if __name__ == '__main__':
    main()
