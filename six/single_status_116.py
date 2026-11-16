#!/usr/bin/env python3
"""
Status of each event in each single comparison of two groups of samples (dataset pairs, fibroblasts, sex / age /
tissue controls), six-dataset version:
  differentially spliced (DS): LeafCutter p.adjust <= 0.05 and max|dPSI| >= 0.1 (as in effect_size_analysis.py);
  splicing conserved: LeafCutter P (not adjusted) > 0.1 and max|dPSI| <= 0.05;
  inconclusive: all other testable events.
max|dPSI| = largest |difference of the group mean PSI| over the junctions of the event (effect_size_analysis.py).
The 'conserved' class is written as 'equivalent' in the status columns (status_m0.025 / 0.05 / 0.1 all carry the new
status), so that the downstream scripts are unchanged; the outputs relabel it 'splicing conserved'.
Writes events/id50_6ds_anova/{effects, effects_confounders}/ (real directories: event_effects.tsv rewritten, the other
files linked to events/id50_6ds/), then summarize_comparisons_116.py.
"""
import os
import subprocess
import sys

import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/'
SRC, DST = B + 'id50_6ds/', B + 'id50_6ds_anova/'
CONS_P, CONS_DPSI = 0.1, 0.05


def raw_p(sub):
    out = []
    root = SRC + f'leafcutter/{sub}/'
    for comp in os.listdir(root):
        if not os.path.isdir(root + comp):
            continue
        for cell in os.listdir(root + comp):
            f = root + f'{comp}/{cell}/leafcutter_ds_cluster_significance.txt'
            if os.path.exists(f):
                d = pd.read_csv(f, sep='\t', usecols=['cluster', 'p'])
                out.append(d.assign(cluster=d['cluster'].str.split(':').str[1], comparison=comp, cell=cell))
    return pd.concat(out).rename(columns={'p': 'p_raw'})


def main():
    for name, sub in (('effects', 'original'), ('effects_confounders', 'confounders')):
        dst = DST + name + '/'
        if os.path.islink(DST + name):
            os.remove(DST + name)
        os.makedirs(dst, exist_ok=True)
        for f in os.listdir(SRC + name):
            if f not in ('event_effects.tsv', 'comparison_summary.tsv', 'comparison_type_summary.tsv') and \
                    not os.path.lexists(dst + f):
                os.symlink(SRC + name + '/' + f, dst + f)
        ev = pd.read_csv(SRC + name + '/event_effects.tsv', sep='\t')
        ev = ev.merge(raw_p(sub), on=['cluster', 'comparison', 'cell'], how='left')
        ds = ev['status_m0.1'] == 'DS'
        cons = ev['testable'] & ~ds & (ev['p_raw'] > CONS_P) & (ev['max_abs_dpsi'] <= CONS_DPSI)
        new = pd.Series('inconclusive', index=ev.index).mask(ds, 'DS').mask(cons, 'equivalent') \
            .mask(ev['status_m0.1'] == 'not_testable', 'not_testable')
        for m in ('0.025', '0.05', '0.1'):
            ev[f'status_ci_m{m}'] = ev[f'status_m{m}']          # the CI-based equivalence status, kept for reference
            ev[f'status_m{m}'] = new
        ev.to_csv(dst + 'event_effects.tsv', sep='\t', index=False)
        t = ev[ev['status_m0.1'] != 'not_testable']
        print(name, t['status_m0.1'].value_counts().to_dict())
        if name == 'effects':
            subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                         'summarize_comparisons_116.py'), dst], check=True,
                           stdout=subprocess.DEVNULL)


if __name__ == '__main__':
    main()
