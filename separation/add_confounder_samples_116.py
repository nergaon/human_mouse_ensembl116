#!/usr/bin/env python3
"""
Confounder control comparisons (reviewers: age, sex, tissue): append the samples to the sample sheet (new rows only,
so existing row numbers / SGE task ids do not change) and write the comparison spec used by confounders_116.py.

Sources (ensembl116/metadata/confounder_comparisons_candidates.tsv, 'best' rows, from a GEO/ENA/ENCODE search):
  human_sex     GSE107011 (Monaco 2019): 2 male vs 2 female donors, sorted B naive, CD4/CD8 naive, NK, classical monocytes
  human_tissue  GSE133383: CD56dim NK cells, blood vs spleen of the same 3 organ donors
  human_age     GSE143213: naive CD4 T cells, young (23-35 y) vs old (65-78 y); G486 dropped (4.8M read pairs)
  mouse_sex     ENCODE spleen B cells (B6CASTF1/J), 8 female vs 8 male
  mouse_tissue  GSE233347: Ly49H+ NK cells, spleen vs liver of the same mice
  mouse_age     GSE184496: naive CD4 T cells, 7-11 weeks vs 83-86 weeks (3 male + 3 female per age)
  and, already in the data:
  mouse_tissue_m1  GSE116177 (M1): bone marrow vs blood monocytes (2 vs 2)
  (GSE124829, ImmGen 2x25, was processed but dropped: too few testable events)
"""
import re

import pandas as pd

MD = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/metadata/'
SHEET = MD + 'sample_sheet_new_datasets.tsv'
CAND = MD + 'confounder_comparisons_candidates.tsv'
SPEC = MD + 'confounder_comparisons.tsv'
STUDY = {'mouse_sex': 'ENCODEB'}  # ENCODE B cells kept apart from the ENCODE4 CD4 T cells of M2
SPECIES = {'human': 'human', 'mouse': 'mouse'}


def tag(r):
    """Short, unique sample tag from the sample title."""
    t = r['sample_title']
    if r['comparison'] == 'mouse_sex':
        m = re.match(r'(ENCSR\w+) rep(\d)', t)
        return f"{r['group'][0].upper()}_{m.group(1)}_{m.group(2)}"
    t = re.sub(r'^NK cells, ', '', t)
    t = re.sub(r'[^A-Za-z0-9]+', '_', t).strip('_')
    return f"{r['group']}_{t}"


def main():
    sheet = pd.read_csv(SHEET, sep='\t', dtype=str, keep_default_na=False)
    cand = pd.read_csv(CAND, sep='\t', dtype=str, keep_default_na=False)
    best = cand[(cand['choice'] == 'best') & ~cand['sample_title'].str.contains('G486')].copy()
    best['study'] = best.apply(lambda r: STUDY.get(r['comparison'], r['accession']), axis=1)
    best['sample'] = best.apply(lambda r: f"{r['cell_type']}_{tag(r)}_{r['study']}", axis=1)
    assert best['sample'].is_unique
    rows = []
    for _, r in best.iterrows():
        species = r['comparison'].split('_')[0]
        runs = r['runs']
        if runs.startswith('ENCFF'):
            fq1 = ';'.join(f'https://www.encodeproject.org/files/{f}/@@download/{f}.fastq.gz' for f in runs.split(','))
        else:
            fq1 = ''
        rows.append(dict(dataset='HC' if species == 'human' else 'MC', study=r['study'], species=species,
                         cell=r['cell_type'], sample=r['sample'], runs=runs, layout=r['layout'],
                         read_count=r['read_count'], fastq_1=fq1, fastq_2=''))
    new = pd.DataFrame(rows)
    new = new[~new['sample'].isin(sheet['sample'])]
    first = len(sheet) + 1
    pd.concat([sheet, new], ignore_index=True).to_csv(SHEET, sep='\t', index=False)
    print(f'appended {len(new)} samples as rows {first}-{first + len(new) - 1}')

    # comparison spec: sample = column name in the counts table ('<sample>_<study>')
    spec = [dict(comparison=r['comparison'], factor=r['comparison'].split('_')[1], species=r['comparison'].split('_')[0],
                 cell=r['cell_type'], group=r['group'], study=r['study'], sample=f"{r['sample']}_{r['study']}")
            for _, r in best.iterrows()]
    for tissue, prefix in (('bone_marrow', 'MonoBM'), ('blood', 'MonoPB')):
        for k in (1, 2):
            spec.append(dict(comparison='mouse_tissue_m1', factor='tissue', species='mouse', cell='Mono',
                             group=tissue, study='GSE116177', sample=f'{prefix}-{k}_GSE116177'))
    spec = pd.DataFrame(spec)
    spec.to_csv(SPEC, sep='\t', index=False)
    print(spec.groupby(['comparison', 'cell', 'group']).size().unstack(fill_value=0).to_string())


if __name__ == '__main__':
    main()
