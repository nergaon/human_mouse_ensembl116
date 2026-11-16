#!/usr/bin/env python3
"""
Orthologous alternative splicing events for the 3 human x 3 mouse design (each dataset >= 3 replicates per cell type).

Human: H1 GSE60424, H2 GSE115736, H3 GSE118165 (as in the main analysis).
Mouse (one study per cell type; dataset_plan_3x3.xlsx):
          B           CD4 T               CD8 T      NK         Monocytes
   M1     GSE287527   GSE215775           GSE313371  GSE261448  GSE261448
   M2     GSE307667   GSE184496 (young)   GSE263132  GSE236819  GSE254202
   M3     GSE231772   GSE134322           GSE284441  GSE180053  PRJNA797236 (blood)
Steps: junction table (leafcutter_JSR_table.py) and mapping to human coordinates (filter_ortholog_junctions.py, id50
orthology) for the new studies; samples not in the design are dropped (GSE184496 old mice, PRJNA797236 spleen);
events with shared_junctions_116.py (same rules as the main analysis: present in >= 2 human and >= 2 mouse
datasets, >= 10 reads in >= 2 samples in >= 2 datasets); sample design (dataset, cell) for the DSA analysis.
Output: ensembl116/events/id50_6ds/ (AS_clusters_value.txt, sample_design.tsv)
"""
import os
import re
import subprocess

import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
E = B + 'ensembl116/'
M3 = E + 'm3/'
J = E + 'junctions/id50/'
OUT = E + 'events/id50_6ds/'
PY = '/gpfs0/tals/projects/software/anaconda3/envs/spyder-env/bin/python'   # networkx
S = B + 'scripts/ensembl116/'
NEW = ['GSE307667', 'GSE231772', 'GSE215775', 'GSE134322', 'GSE263132', 'GSE284441', 'GSE236819', 'GSE180053',
       'GSE254202']
MOUSE = {  # label -> (study, table, cell, sample filter)
    'M1': [('GSE287527', J + 'GSE287527_orthologs_junctions.tsv', 'NveB', None),
           ('GSE215775', M3 + 'junctions/GSE215775_orthologs_junctions.tsv', 'CD4T', None),
           ('GSE313371', J + 'GSE313371_orthologs_junctions.tsv', 'CD8T', None),
           ('GSE261448', J + 'GSE261448_orthologs_junctions.tsv', None, None)],
    'M2': [('GSE307667', M3 + 'junctions/GSE307667_orthologs_junctions.tsv', 'NveB', None),
           ('GSE184496', J + 'GSE184496_orthologs_junctions.tsv', 'CD4T', 'young'),
           ('GSE263132', M3 + 'junctions/GSE263132_orthologs_junctions.tsv', 'CD8T', None),
           ('GSE236819', M3 + 'junctions/GSE236819_orthologs_junctions.tsv', 'NK', None),
           ('GSE254202', M3 + 'junctions/GSE254202_orthologs_junctions.tsv', 'Mono', None)],
    'M3': [('GSE231772', M3 + 'junctions/GSE231772_orthologs_junctions.tsv', 'NveB', None),
           ('GSE134322', M3 + 'junctions/GSE134322_orthologs_junctions.tsv', 'CD4T', None),
           ('GSE284441', M3 + 'junctions/GSE284441_orthologs_junctions.tsv', 'CD8T', None),
           ('GSE180053', M3 + 'junctions/GSE180053_orthologs_junctions.tsv', 'NK', None),
           ('PRJNA797236', J + 'PRJNA797236_orthologs_junctions.tsv', 'Mono', 'blood')],
}
HUMAN = {'H1': J + 'GSE60424_orthologs_junctions.tsv', 'H2': J + 'GSE115736_orthologs_junctions.tsv',
         'H3': J + 'GSE118165_orthologs_junctions.tsv'}


def run(cmd, cwd=None):
    print('+', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=cwd, stdout=subprocess.DEVNULL)


def main():
    os.makedirs(M3 + 'junctions', exist_ok=True)
    os.makedirs(OUT + 'tables', exist_ok=True)
    for st in NEW:   # junction table and mapping of the new studies
        d = M3 + f'{st}/leafcutter_0.2.9/'
        juncs = sorted(f for f in os.listdir(d) if f.endswith('.junc'))
        open(d + 'juncfiles.txt', 'w').write('\n'.join(juncs) + '\n')
        run(['python3', S + 'leafcutter_JSR_table.py', '-j', 'juncfiles.txt', '-o', f'{st}_JSR', '-l', '100000',
             '--nochromcheck'], cwd=d)
        run(['python3', S + 'filter_ortholog_junctions.py', '--points', E + 'orthology/unique_points_116_id50.txt',
             '--mode', 'remap', '--source-column', 'position_m', '--target-column', 'position_h',
             '--inputs', d + f'{st}_JSR_junction_counts_no_strand.tsv',
             '--outputs', M3 + f'junctions/{st}_orthologs_junctions.tsv'])
    args = []
    for lab, path in HUMAN.items():
        args += ['--table', f'{lab}:human:{path}']
    for lab, items in MOUSE.items():
        for st, path, cell, filt in items:
            t = pd.read_csv(path, sep='\t', index_col=0)
            keep = [c for c in t.columns if c == 'original_junction_id' or
                    ((cell is None or c.startswith(cell + '_')) and (filt is None or f'_{filt}_' in c))]
            sub = OUT + f'tables/{lab}_{st}.tsv'
            t[keep].to_csv(sub, sep='\t')
            args += ['--table', f'{lab}:mouse:{sub}']
    run([PY, S + 'shared_junctions_116.py'] + args + ['--out-dir', OUT])
    # sample design: human samples as in the main analysis; mouse samples by label and cell prefix
    cols = pd.read_csv(OUT + 'AS_clusters_value.txt', sep='\t', index_col=0, nrows=1).columns
    old = pd.read_csv(E + 'events/id50_6ds/leafcutter/sample_design.tsv', sep='\t')
    old = old[old['dataset'].isin(['H1', 'H2', 'H3'])].set_index('sample')
    rows = []
    for c in cols:
        if c in old.index:
            rows.append((c, old.loc[c, 'dataset'], old.loc[c, 'cell']))
            continue
        for lab, items in MOUSE.items():
            if not c.endswith('_' + lab):   # columns are '<table>.tsv_<sample>_<label>'
                continue
            for st, path, cell, filt in items:
                if st in c:
                    m = re.search(r'(?:^|_)(NveB|CD4T|CD8T|NK|Mono)_', c)
                    if m and (cell is None or m.group(1) == cell) and (filt is None or f'_{filt}_' in c):
                        rows.append((c, lab, m.group(1)))
    design = pd.DataFrame(rows, columns=['sample', 'dataset', 'cell']).drop_duplicates('sample')
    design.to_csv(OUT + 'sample_design.tsv', sep='\t', index=False)
    print(design.groupby(['cell', 'dataset']).size().unstack(fill_value=0))


if __name__ == '__main__':
    main()
