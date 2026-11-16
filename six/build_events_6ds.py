#!/usr/bin/env python3
"""
Orthologous alternative splicing events of the 3 human x 3 mouse design (each dataset >= 3 replicates per cell type),
for one orthology variant, and the LeafCutter design of every dataset pair (id50 only).

Human: H1 GSE60424 (healthy controls), H2 GSE115736, H3 GSE118165 (unstimulated).
Mouse (one study per cell type; ensembl116/m3/dataset_plan_3x3.xlsx):
          B           CD4 T               CD8 T      NK         Monocytes
   M1     GSE287527   GSE215775           GSE313371  GSE261448  GSE261448
   M2     GSE307667   GSE184496 (young)   GSE263132  GSE236819  GSE254202
   M3     GSE231772   GSE134322           GSE284441  GSE180053  PRJNA797236 (blood)
Fibroblasts (counts only, not used to define the events): E-MTAB-5919 human (FH) and mouse (FM).
Steps: mouse junctions remapped to human coordinates with the variant's orthologous boundaries (studies not yet mapped
for the variant); per study, only the samples of the design; events with shared_junctions_116.py (present in >= 2
human and >= 2 mouse datasets, >= 10 reads in >= 2 samples in >= 2 datasets); sample design (dataset, cell);
LeafCutter groups for the 15 dataset pairs per cell type (9 human-mouse, 3 human-human, 3 mouse-mouse) and FH-FM.
Output: ensembl116/events/<variant>_6ds/ (AS_clusters_value.txt, sample_design.tsv, leafcutter/)
Usage: build_events_6ds.py <variant>   (id50, id70, id80, id50_rbh)
"""
import itertools
import os
import re
import subprocess
import sys

import pandas as pd

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
E = B + 'ensembl116/'
S = B + 'scripts/ensembl116/'
PY = '/gpfs0/tals/projects/software/anaconda3/envs/spyder-env/bin/python'   # networkx
CELLS = ['NveB', 'CD4T', 'CD8T', 'NK', 'Mono']
HUMAN = {'H1': 'GSE60424', 'H2': 'GSE115736', 'H3': 'GSE118165'}
MOUSE = {  # label -> [(study, cell or None = all cells, sample filter)]
    'M1': [('GSE287527', 'NveB', None), ('GSE215775', 'CD4T', None), ('GSE313371', 'CD8T', None),
           ('GSE261448', None, None)],
    'M2': [('GSE307667', 'NveB', None), ('GSE184496', 'CD4T', 'young'), ('GSE263132', 'CD8T', None),
           ('GSE236819', 'NK', None), ('GSE254202', 'Mono', None)],
    'M3': [('GSE231772', 'NveB', None), ('GSE134322', 'CD4T', None), ('GSE284441', 'CD8T', None),
           ('GSE180053', 'NK', None), ('PRJNA797236', 'Mono', 'blood')],
}
NEW = ['GSE307667', 'GSE231772', 'GSE215775', 'GSE134322', 'GSE263132', 'GSE284441', 'GSE236819', 'GSE180053',
       'GSE254202']
JSR = {st: E + f'm3/{st}/leafcutter_0.2.9/{st}_JSR_junction_counts_no_strand.tsv' for st in NEW}
JSR['PRJNA797236'] = B + 'PRJNA797236/leafcutter_0.2.9/PRJNA797236_JSR_junction_counts_no_strand.tsv'
FIBRO = {'FH': ('human', 'EMTAB5919H', r'^HS\d+_UNST'), 'FM': ('mouse', 'EMTAB5919M', r'^MM\d+_UNST')}
H_CELL = {'GSE60424': [(r'^Bcells_HealthyControl', 'NveB'), (r'^CD4_HealthyControl', 'CD4T'),
                       (r'^CD8_HealthyControl', 'CD8T'), (r'^NK_HealthyControl', 'NK'),
                       (r'^Monocytes_HealthyControl', 'Mono')],
          'GSE115736': [(rf'^{c}_\d+_', c) for c in CELLS],
          'GSE118165': [(rf'^{c}_', c) for c in CELLS]}


def run(cmd, cwd=None):
    print('+', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=cwd, stdout=subprocess.DEVNULL)


def junction_table(study, v):
    """Orthologous junction table of a study for variant v (mapped here if needed)."""
    p = E + f'junctions/{v}/{study}_orthologs_junctions.tsv'
    if os.path.exists(p):
        return p
    d = E + f'm3/junctions/{v}/'
    p = d + f'{study}_orthologs_junctions.tsv'
    if v == 'id50' and os.path.exists(E + f'm3/junctions/{study}_orthologs_junctions.tsv'):
        return E + f'm3/junctions/{study}_orthologs_junctions.tsv'
    if not os.path.exists(p):
        os.makedirs(d, exist_ok=True)
        if not os.path.exists(JSR[study]):   # JSR table of a new study
            jd = os.path.dirname(JSR[study]) + '/'
            open(jd + 'juncfiles.txt', 'w').write('\n'.join(sorted(f for f in os.listdir(jd)
                                                                   if f.endswith('.junc'))) + '\n')
            run(['python3', S + 'leafcutter_JSR_table.py', '-j', 'juncfiles.txt', '-o', f'{study}_JSR', '-l',
                 '100000', '--nochromcheck'], cwd=jd)
        run(['python3', S + 'filter_ortholog_junctions.py', '--points', E + f'orthology/unique_points_116_{v}.txt',
             '--mode', 'remap', '--source-column', 'position_m', '--target-column', 'position_h',
             '--inputs', JSR[study], '--outputs', p])
    return p


def main():
    v = sys.argv[1]
    out = E + f'events/{v}_6ds/'
    os.makedirs(out + 'tables', exist_ok=True)
    args = []
    for lab, st in HUMAN.items():
        args += ['--table', f'{lab}:human:{junction_table(st, v)}']
    for lab, items in MOUSE.items():
        for st, cell, filt in items:
            t = pd.read_csv(junction_table(st, v), sep='\t', index_col=0)
            keep = [c for c in t.columns if c == 'original_junction_id' or
                    ((cell is None or c.startswith(cell + '_')) and (filt is None or f'_{filt}_' in c))]
            sub = out + f'tables/{st}_orthologs_junctions.tsv'
            t[keep].to_csv(sub, sep='\t')
            args += ['--table', f'{lab}:mouse:{sub}']
    for lab, (sp, st, _) in FIBRO.items():
        args += ['--extra', f'{lab}:{sp}:{junction_table(st, v)}']
    run([PY, S + 'shared_junctions_116.py'] + args + ['--out-dir', out])

    # sample design
    cols = pd.read_csv(out + 'AS_clusters_value.txt', sep='\t', index_col=0, nrows=1).columns
    lab_of = {st: lab for lab, items in MOUSE.items() for st, _, _ in items}
    lab_of.update({st: lab for lab, st in HUMAN.items()})
    rows = []
    for c in cols:
        sample, study = c.rsplit('_', 1)
        if study in HUMAN.values():
            cell = next((cl for pat, cl in H_CELL[study] if re.search(pat, c)), None)
        elif study in lab_of:
            m = re.match(r'(NveB|CD4T|CD8T|NK|Mono)_', sample)
            cell = m.group(1) if m else None
        else:
            lab = next((k for k, (_, s, pat) in FIBRO.items() if s == study and re.search(pat, sample)), None)
            if lab:
                rows.append((c, lab, 'Fibro'))
            continue
        if cell:
            rows.append((c, lab_of[study], cell))
    design = pd.DataFrame(rows, columns=['sample', 'dataset', 'cell'])
    design['matched'] = True
    os.makedirs(out + 'leafcutter', exist_ok=True)
    design.to_csv(out + 'sample_design.tsv', sep='\t', index=False)
    design.to_csv(out + 'leafcutter/sample_design.tsv', sep='\t', index=False)
    print(design.groupby(['cell', 'dataset']).size().unstack(fill_value=0))
    print('not used:', [c for c in cols if c not in set(design['sample'])])
    if v != 'id50':
        return
    hum, mou = list(HUMAN), list(MOUSE)
    pairs = list(itertools.product(hum, mou)) + list(itertools.combinations(hum, 2)) + \
        list(itertools.combinations(mou, 2))
    jobs = []
    for a, b in pairs + [('FH', 'FM')]:
        for cell in CELLS + ['Fibro']:
            ga = design.loc[(design['dataset'] == a) & (design['cell'] == cell), 'sample']
            gb = design.loc[(design['dataset'] == b) & (design['cell'] == cell), 'sample']
            if len(ga) < 2 or len(gb) < 2:
                continue
            path = out + f'leafcutter/original/{a}_{b}/{cell}'
            os.makedirs(path, exist_ok=True)
            pd.concat([pd.DataFrame({'s': ga, 'g': a}), pd.DataFrame({'s': gb, 'g': b})]).to_csv(
                path + '/groups_file.txt', sep='\t', header=False, index=False)
            jobs.append(path)
    open(out + 'leafcutter/jobs.txt', 'w').write('\n'.join(jobs) + '\n')
    print(len(jobs), 'LeafCutter runs')


if __name__ == '__main__':
    main()
