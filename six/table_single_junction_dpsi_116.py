#!/usr/bin/env python3
"""
Candidate table: per immune cell type, the number of events and the number of events in which exactly one junction
(or at least one junction) differs between human and mouse by |dPSI| >= 0.1 or >= 0.2.

dPSI H-M of a junction = mean PSI of H1, H2, H3 - mean PSI of M1, M2 (dataset PSI as in the classification).
Events: those with a PSI in all five datasets (testable in the classification).
Note: in an event with two junctions the two PSIs sum to 1, so both junctions have the same |dPSI|; exactly one
junction above a cutoff is possible only in events with three or more junctions.

Output: ensembl116/revision_6ds/candidate_figures/single_junction_dpsi.(tsv|xlsx)
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import classify_separation_116 as cs  # noqa: E402
import build_supp_tables_116 as bst  # noqa: E402

E = cs.BASE + 'events/id50_6ds/'
OUT = cs.BASE + 'revision_6ds/candidate_figures/'
CELL_NAME = {'CD4T': 'CD4 T', 'CD8T': 'CD8 T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Monocytes'}


def main():
    counts = pd.read_csv(E + 'AS_clusters_value.txt', sep='\t', index_col=0)
    design = pd.read_csv(E + 'leafcutter/sample_design.tsv', sep='\t').set_index('sample')
    rows = []
    for cell in cs.CELLS:
        psi, clu = cs.dataset_psi(counts, design, cell)
        d = (psi[list(cs.HUMAN)].mean(axis=1) - psi[['M1', 'M2']].mean(axis=1)).abs()
        d = pd.DataFrame({'clu': clu, 'd': d.to_numpy()})
        g = d.groupby('clu')['d']
        nj = g.size()
        row = {'Cell type': CELL_NAME[cell], 'Events': len(nj), 'Events with >= 3 junctions': int((nj >= 3).sum())}
        for cut in (0.1, 0.2):
            n = g.apply(lambda v: int((v >= cut).sum()))
            row[f'>= 1 junction with |dPSI| >= {cut}'] = int((n >= 1).sum())
            row[f'Exactly 1 junction with |dPSI| >= {cut}'] = int((n == 1).sum())
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(OUT + 'single_junction_dpsi.tsv', sep='\t', index=False)
    bst.OUT = OUT
    bst.write('single_junction_dpsi', 'Events with one junction that differs between human and mouse',
              'dPSI H-M of a junction = mean PSI of H1, H2, H3 minus mean PSI of M1, M2 (dataset PSI = mean over '
              'samples with >= 20 reads in the event). Events: PSI in all five datasets. In an event with two '
              'junctions both junctions have the same |dPSI|, so exactly one junction above a cutoff requires three or '
              'more junctions.', {'by_cell_type': t})
    print(t.to_string(index=False))


if __name__ == '__main__':
    main()
