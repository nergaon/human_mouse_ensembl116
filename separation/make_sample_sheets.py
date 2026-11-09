#!/usr/bin/env python3
"""
Sample sheets for the datasets added in the revision (FASTQ URLs from ENA / ENCODE).

  H3  GSE118165 (Calderon 2019), human blood, unstimulated (-U) naive B / CD4 / CD8, NK, monocytes
  M2  composite mouse set, one study per cell type:
        NveB  GSE287527 wild-type spleen follicular B cells
        NK    GSE261448 PBS-control spleen NK cells
        Mono  GSE261448 PBS-control bone-marrow Ly6C+ monocytes
        CD4T  ENCODE4 naive CD4 T cells (spleen, B6CASTF1)
        CD8T  GSE313371 uninfected wild-type (B6129) spleen CD8 T cells
              (ENCODE4 has a single naive CD8 replicate; it stays in the sheet but is not used)
  PBMC comparison (separate cell type)
  H4  GSE107011 (Monaco 2019) PBMC, 13 healthy donors
  M3  GSE182858 PBMC of mice co-housed with pet-store mice, baseline (CoH_0_Naive), 13 mice

New rows are appended to an existing sheet so that row numbers (= SGE task ids) never change.
Output: ensembl116/metadata/sample_sheet_new_datasets.tsv with columns
  dataset, study, species, cell, sample, runs, layout, read_count, fastq_1, fastq_2
(fastq_1 / fastq_2 are ';'-separated URL lists; fastq_2 is empty for single-end)
"""
import io
import os
import re

import pandas as pd
import requests

OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/metadata/sample_sheet_new_datasets.tsv'
ENA = ('https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run'
       '&fields=run_accession,sample_title,library_layout,read_count,fastq_ftp&accession=')

# study, bioproject, dataset label, species, title regex -> cell, sample name from regex groups
ENA_SETS = [
    ('GSE118165', 'PRJNA484735', 'H3', 'human',
     r'^(\d+)-(Naive_B|Naive_Teffs|Naive_CD8_T|Mature_NK|Monocytes)-U$',
     {'Naive_B': 'NveB', 'Naive_Teffs': 'CD4T', 'Naive_CD8_T': 'CD8T', 'Mature_NK': 'NK', 'Monocytes': 'Mono'},
     lambda m: (m.group(2), m.group(1))),
    ('GSE287527', 'PRJNA1213147', 'M2', 'mouse', r'^Follicular_B_cells_wt_rep(\d+)$',
     {'B': 'NveB'}, lambda m: ('B', m.group(1))),
    ('GSE261448', 'PRJNA1086964', 'M2', 'mouse', r'^(\d)(NK|Monocytes)-Control$',
     {'NK': 'NK', 'Monocytes': 'Mono'}, lambda m: (m.group(2), m.group(1))),
    ('GSE313371', 'PRJNA1379344', 'M2', 'mouse', r'^B6129, Uninfected, (\d)$',
     {'CD8': 'CD8T'}, lambda m: ('CD8', m.group(1))),
    # confounder controls (reviewers: age, sex): ImmGen sex study, spleen, C57BL/6J; 'F' = female, unlabeled = male
    #   sex arm: 6-week-old female vs male (#1-#3); age arm: 2 months vs 17 (F) / 20 (M) months, one mouse per sex
    ('GSE124829', 'PRJNA513879', 'MC', 'mouse',
     r'^(B|NK|T\.?4|T\.?8)\.(?:(2m|17m|20m)\.(F\.)?|(F\.|M\.))Sp#([123])$',
     {'B': 'NveB', 'NK': 'NK', 'T.4': 'CD4T', 'T4': 'CD4T', 'T.8': 'CD8T', 'T8': 'CD8T'},
     lambda m: (m.group(1), f"{m.group(2) or '6w'}_{'F' if (m.group(3) or m.group(4)) == 'F.' else 'M'}_{m.group(5)}")),
    ('GSE107011', 'PRJNA418779', 'H4', 'human', r'^(\w+)_PBMC_rep\d+$',
     {'PBMC': 'PBMC'}, lambda m: ('PBMC', m.group(1))),
    ('GSE182858', 'PRJNA758072', 'M3', 'mouse', r'^(\d+)_CoH_0_Naive_S\d+ \[PBMC vaccine transcriptome\]$',
     {'PBMC': 'PBMC'}, lambda m: ('PBMC', m.group(1))),
]
ENCODE_SETS = {'CD4T': ['ENCSR968ZOA', 'ENCSR835KFT'], 'CD8T': ['ENCSR769LJY']}


def ena_rows():
    rows = []
    for study, project, dataset, species, pattern, cells, key in ENA_SETS:
        runs = pd.read_csv(io.StringIO(requests.get(ENA + project, timeout=120).text), sep='\t')
        for title, g in runs.groupby('sample_title'):
            m = re.match(pattern, title)
            if not m:
                continue
            label, rep = key(m)
            cell = cells[label]
            # ENA FASTQ links are informational: process_sample.sh fetches SRR runs from the NCBI AWS mirror
            urls = [u.split(';') for u in g['fastq_ftp'].fillna('') if u]
            layout = g['library_layout'].iloc[0]
            r1 = ['https://' + next(x for x in u if x.endswith('_1.fastq.gz') or len(u) == 1) for u in urls]
            r2 = ['https://' + next(x for x in u if x.endswith('_2.fastq.gz')) for u in urls] if layout == 'PAIRED' else []
            rows.append(dict(dataset=dataset, study=study, species=species, cell=cell,
                             sample=f'{cell}_{rep}_{study}', runs=','.join(g['run_accession']), layout=layout,
                             read_count=int(g['read_count'].sum()), fastq_1=';'.join(r1), fastq_2=';'.join(r2)))
    return rows


def encode_rows():
    rows = []
    for cell, experiments in ENCODE_SETS.items():
        for k, exp in enumerate(experiments, 1):
            url = (f'https://www.encodeproject.org/search/?type=File&dataset=/experiments/{exp}/'
                   '&file_format=fastq&status=released&format=json&limit=100')
            files = requests.get(url, headers={'Accept': 'application/json'}, timeout=120).json()['@graph']
            single = [f for f in files if f.get('run_type') == 'single-ended']
            rows.append(dict(dataset='M2', study='ENCODE4', species='mouse', cell=cell, sample=f'{cell}_{k}_ENCODE4',
                             runs=','.join(exp + ':' + f['accession'] for f in single), layout='SINGLE',
                             read_count=sum(f.get('read_count') or 0 for f in single),
                             fastq_1=';'.join('https://www.encodeproject.org' + f['href'] for f in single), fastq_2=''))
    return rows


def main():
    sheet = pd.DataFrame(ena_rows() + encode_rows()).sort_values(['dataset', 'cell', 'sample'])
    if os.path.exists(OUT):  # keep existing rows and their order; append new samples
        old = pd.read_csv(OUT, sep='\t', dtype=str, keep_default_na=False)
        sheet = pd.concat([old, sheet[~sheet['sample'].isin(old['sample'])].astype(str)], ignore_index=True)
    sheet.to_csv(OUT, sep='\t', index=False)
    print(len(sheet), 'rows')
    print(sheet.groupby(['dataset', 'study', 'cell']).size())


if __name__ == '__main__':
    main()
