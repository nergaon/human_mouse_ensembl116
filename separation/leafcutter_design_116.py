#!/usr/bin/env python3
"""
LeafCutter comparison design: groups files for every dataset pair x cell type.

Datasets
  H1 GSE60424 (healthy controls)   H2 GSE115736   H3 GSE118165 (unstimulated)
  M1 GSE116177                      M2 composite (GSE287527 B, GSE261448 NK/Mono, ENCODE4 CD4, GSE313371 CD8)
  Fibroblasts only: FH / FM E-MTAB-5919 (human / mouse dermal fibroblasts)
Comparisons
  cross-species: every human x mouse pair (H2-M1 is the pair used in the manuscript)
  within-species: every human pair, M1-M2
Sample sets
  original  M1 groups as in the manuscript (CD4T incl. effector/memory, CD8T incl. memory,
            Mono incl. inflammatory and bone-marrow monocytes)
  matched   M1 restricted to naive CD4/CD8 T and blood monocytes (reviewer 4)
Cell types without >= 2 samples in both groups are skipped.

Output: <out-dir>/<set>/<A>_<B>/<cell>/groups_file.txt and <out-dir>/jobs.txt (one directory per line)
Usage: leafcutter_design_116.py --counts AS_clusters_value.txt --out-dir DIR
"""
import argparse
import itertools
import os
import re

import pandas as pd

CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono', 'Fibro']
HUMAN, MOUSE = ['H1', 'H2', 'H3', 'FH'], ['M1', 'M2', 'FM']

# study -> list of (regex on the sample name, cell, matched-set member)
RULES = {
    'GSE60424': [(r'^Bcells_HealthyControl', 'NveB', True), (r'^CD4_HealthyControl', 'CD4T', True),
                 (r'^CD8_HealthyControl', 'CD8T', True), (r'^NK_HealthyControl', 'NK', True),
                 (r'^Monocytes_HealthyControl', 'Mono', True)],
    'GSE115736': [(rf'^{c}_\d+$', c, True) for c in CELLS],
    'GSE116177': [(r'^BCell-', 'NveB', True), (r'^NveCD4T-', 'CD4T', True), (r'^(CD4T|EffCD4T|MemCD4T)-', 'CD4T', False),
                  (r'^NveCd8T-', 'CD8T', True), (r'^(CD8T|MemCD8T)-', 'CD8T', False), (r'^NK-', 'NK', True),
                  (r'^MonoPB-', 'Mono', True), (r'^(InfMono|MonoBM)-', 'Mono', False)],
}
for study in ('GSE118165', 'GSE287527', 'GSE261448', 'GSE313371'):
    RULES[study] = [(rf'^{c}_', c, True) for c in CELLS]
RULES['ENCODE4'] = [(r'^CD4T_', 'CD4T', True)]  # the single ENCODE4 CD8 replicate is not used
RULES['EMTAB5919H'] = [(r'^HS\d+_UNST$', 'Fibro', True)]
RULES['EMTAB5919M'] = [(r'^MM\d+_UNST$', 'Fibro', True)]
DATASET = {'GSE60424': 'H1', 'GSE115736': 'H2', 'GSE118165': 'H3', 'GSE116177': 'M1',
           'GSE287527': 'M2', 'GSE261448': 'M2', 'ENCODE4': 'M2', 'GSE313371': 'M2',
           'EMTAB5919H': 'FH', 'EMTAB5919M': 'FM'}


def assign(column):
    """'<sample>_<study>' column of the LeafCutter count table -> (dataset, cell, matched) or None."""
    sample, study = column.rsplit('_', 1)
    for pattern, cell, matched in RULES.get(study, []):
        if re.search(pattern, sample):
            return DATASET[study], cell, matched
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    columns = pd.read_csv(args.counts, sep='\t', index_col=0, nrows=1).columns
    design = pd.DataFrame([(c, *a) for c in columns if (a := assign(c))],
                          columns=['sample', 'dataset', 'cell', 'matched'])
    unassigned = [c for c in columns if assign(c) is None]
    print(f'{len(design)} samples assigned; not used: {unassigned}')
    print(design.groupby(['dataset', 'cell', 'matched']).size().unstack(fill_value=0))
    design.to_csv(os.path.join(args.out_dir, 'sample_design.tsv'), sep='\t', index=False)

    pairs = list(itertools.product(HUMAN, MOUSE)) + list(itertools.combinations(HUMAN, 2)) + [('M1', 'M2')]
    # pairs or cell types without samples are skipped below
    jobs = []
    for sample_set in ('original', 'matched'):
        d = design if sample_set == 'original' else design[design['matched']]
        for a, b in pairs:
            if sample_set == 'matched' and 'M1' not in (a, b):
                continue  # only M1 groups change
            for cell in CELLS:
                ga = d[(d['dataset'] == a) & (d['cell'] == cell)]['sample']
                gb = d[(d['dataset'] == b) & (d['cell'] == cell)]['sample']
                if len(ga) < 2 or len(gb) < 2:
                    if len(ga) or len(gb):
                        print(f'skip {sample_set} {a}-{b} {cell}: {len(ga)} vs {len(gb)} samples')
                    continue
                path = os.path.join(args.out_dir, sample_set, f'{a}_{b}', cell)
                os.makedirs(path, exist_ok=True)
                pd.concat([pd.DataFrame({'s': ga, 'g': a}), pd.DataFrame({'s': gb, 'g': b})]) \
                    .to_csv(os.path.join(path, 'groups_file.txt'), sep='\t', header=False, index=False)
                jobs.append(path)
    with open(os.path.join(args.out_dir, 'jobs.txt'), 'w') as f:
        f.write('\n'.join(jobs) + '\n')
    print(len(jobs), 'LeafCutter runs')


if __name__ == '__main__':
    main()
