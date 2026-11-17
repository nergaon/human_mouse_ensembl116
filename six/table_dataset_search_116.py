#!/usr/bin/env python3
"""
Supplementary Table S3: the datasets of this study and the candidate datasets found in the searches for additional
public datasets (October 2026), one row per accession, with whether each was used in the main human-mouse analysis
and in the sex, age and tissue control analysis, and why not.

Sources: ensembl116/metadata/candidate_datasets.tsv (immune datasets of the five cell types),
ensembl116/metadata/confounder_comparisons_candidates.tsv (control comparisons), ensembl116/metadata/m3_candidates.tsv
and ensembl116/m3/dataset_plan_3x3.xlsx (per-cell-type search for the mouse datasets of the 3 human x 3 mouse design),
Table S2 (datasets used), and the ENA, GEO and ENCODE records.

Output: revision_6ds/supplementary/Table_S3.xlsx
"""
import re

import openpyxl
import pandas as pd
from openpyxl.styles import Font

import build_supp_tables_116 as bst

MAIN, CTRL = 'Main analysis (human vs mouse)', 'Control analysis (sex, age, tissue)'
COLS = ['Species', 'Accession', 'Study / source', 'Cell types', 'Tissue', 'Biological replicates per cell type',
        'Read layout and length', 'Reads per sample', MAIN, 'Reason not used in the main analysis', CTRL,
        'Reason not used in the control analysis']

H, M = 'Human', 'Mouse'
SAME = 'Same as for the main analysis.'
NOGROUPS = 'Sex and age not given; one tissue.'
ROWS = [
    # ---------------- human
    (H, 'GSE60424 (H1)', 'Linsley et al. 2014; sorted whole-blood immune cells (healthy controls)',
     'B, CD4 T, CD8 T, NK, monocytes', 'Peripheral blood', '4 donors', 'Paired, 2x50', '9-29M pairs',
     'Yes: H1', '', 'No', 'All donors female, aged 22-32 years; one tissue.'),
    (H, 'GSE115736 (H2)', 'Choi et al. 2019; Haemopedia human', 'B, CD4 T, CD8 T, NK, monocytes',
     'Peripheral blood', '5 donors', 'Single, 1x100', '16-22M', 'Yes: H2', '', 'No', NOGROUPS),
    (H, 'GSE118165 (H3)', 'Calderon et al. 2019; sorted blood cell types, resting and stimulated',
     'Naive B, naive CD4 T, naive CD8 T, mature NK, monocytes', 'Peripheral blood', '4 donors',
     'Paired, 2x75-2x100', '8-71M pairs', 'Yes: H3 (unstimulated samples only)', '', 'No', NOGROUPS),
    (H, 'E-MTAB-5919 (FH)', 'Human and mouse dermal fibroblasts', 'Fibroblasts', 'Dermis', '6 donors',
     'Paired, 2x125', '25-36M pairs', 'Yes: fibroblast comparison (FH)', '', 'No',
     'All donors female, aged 40-49 years; one tissue.'),
    (H, 'GSE107011', 'Monaco et al. 2019; 29 sorted PBMC cell types', 'Naive B, naive CD4 T, naive CD8 T, NK, '
     'classical monocytes', 'PBMC', '4 donors (2 male, 2 female)', 'Paired, 2x50', '12-24M pairs', 'No',
     'Used as the sex control (Supplementary Note 1); three human datasets per cell type were sufficient.', 'Yes: human sex (2 male vs 2 female donors, five cell types)', ''),
    (H, 'GSE143213', 'Sorted blood naive CD4 T cells (PMID 32041953)', 'Naive CD4 T', 'Peripheral blood',
     '4-5 per age group', 'Paired, 2x50', '-', 'No', 'Profiles only one of the five cell types.',
     'Yes: human age (young, 23-35 years, vs old, 65-78 years)', ''),
    (H, 'GSE133383', 'Dogra et al. 2020; NK cells of organ donors', 'CD56dim CD16+ NK', 'Blood, spleen, bone marrow',
     '3 donors', 'Paired, 2x150', '36-50M pairs', 'No',
     'Profiles only one of the five cell types; deceased organ donors.',
     'Yes: human tissue (blood vs spleen of the same donors)', ''),
    (H, 'GSE282496', 'Sorted naive CD4 T cells, unstimulated time point', 'Naive CD4 T', 'Cord blood, adult blood',
     '-', 'Paired, 2x150', '-', 'No', 'Profiles only one of the five cell types.', 'No',
     'Neonatal vs adult, not young vs old adults; adult ages not given.'),
    (H, 'GSE94964', 'Kumar et al. 2017; organ and blood donors', 'CD69- memory CD8 T', 'Spleen, blood', '3 per tissue',
     'Single, 1x101', '17-32M', 'No', 'Memory, not naive, CD8 T cells only.', 'No',
     'Tissue fully confounded with donor (spleen from organ donors, blood from other, living donors).'),
    (H, 'GSE74246', 'Corces et al. 2016; hematopoietic cell types and AML', 'B, CD4 T, CD8 T, NK, monocytes',
     'Peripheral blood', '4 donors, of which 2-3 per cell type with 2x76 reads',
     'Paired, 2x76 (5 donors) and 2x36 (2 donors)', '14-66M pairs', 'No',
     'Only 2-3 samples per cell type with reads longer than 36 nt; total (not naive) B and T cells; donors not '
     'matched between cell types; SPIA amplification.', 'No', 'Sex and age not given; one tissue.'),
    (H, 'GSE199324', 'White et al. 2024; DREAM deconvolution challenge', 'Naive B, naive CD4 T, naive CD8 T, NK, '
     'monocytes (+ other types and admixtures)', 'PBMC (commercial vendors)', '1-2 per cell type',
     'Paired, 2x100', 'Deep', 'No', 'Only one or two replicates per cell type (one for naive B and naive CD8 T); '
     'donors not matched.', 'No', 'Sex and age not given; one or two replicates per cell type.'),
    (H, 'GSE64655', 'Hoek et al. 2015; influenza vaccination time course', 'B, total T (not CD4/CD8), NK, monocytes',
     'Peripheral blood', '2 donors (day 0)', 'Paired, 2x50', '19-32M pairs', 'No',
     'Two donors only; T cells not separated into CD4 and CD8; short reads; other time points post-vaccination.',
     'No', SAME),
    (H, 'GSE149050', 'Andreoletti et al. 2021; SLE vs healthy controls', 'T, B, NK and other PBMC-sorted cells',
     'PBMC', '~12 healthy controls', 'Single, 1x50', '~16M', 'No', 'Single-end 50 nt reads, shallow.', 'No', SAME),
    (H, 'DICE (dbGaP phs001703)', 'Schmiedel et al. 2018; 15 sorted immune cell types', 'Naive B, naive CD4 T, '
     'naive CD8 T, NK, classical monocytes', 'Peripheral blood', '91 donors', 'Not public', '-', 'No',
     'Raw reads under controlled access (dbGaP); only gene-level expression is public, no junction reads.',
     'No', SAME),
    (H, 'BLUEPRINT (EGA) / ImmuNexUT (NBDC)', 'Sorted immune cell types', 'Monocytes, naive CD4 T and others',
     'Peripheral blood', 'Many', 'Paired, 2x100 (BLUEPRINT)', '-', 'No', 'Raw reads under controlled access.',
     'No', SAME),
    # ---------------- mouse
    (M, 'E-MTAB-5919 (FM)', 'Human and mouse dermal fibroblasts', 'Fibroblasts', 'Dermis', '3 mice',
     'Paired, 2x125', '25-36M pairs', 'Yes: fibroblast comparison (FM)', '', 'No',
     'All mice female, 8 weeks; one tissue.'),
    # M1, M2 and M3: three mouse datasets per cell type (B, CD4 T, CD8 T, NK cells, monocytes)
    (M, 'GSE287527', 'Mek1/Mek2 deletion in spleen B-cell subsets (wild-type samples)', 'Follicular B',
     'Spleen', '4', 'Paired, 2x101', '56-72M pairs', 'Yes: M1 B cells (wild type)', '', 'No',
     'Sex not given; one age group (3-5 months); one tissue.'),
    (M, 'GSE215775', 'Naive CD4 T cells of spleen and mesenteric lymph nodes (wild type)', 'Naive CD4 T',
     'Spleen (mesenteric lymph nodes not used)', '3 per tissue', 'Paired, 2x151', '24-30M pairs',
     'Yes: M1 CD4 T cells (spleen samples)', '', 'No',
     'Sex and age not given.'),
    (M, 'GSE313371', 'Uninfected wild-type B6129 mice', 'CD8 T', 'Spleen', '3', 'Paired, 2x150', '16-19M pairs',
     'Yes: M1 CD8 T cells (uninfected wild type)', '', 'No', NOGROUPS),
    (M, 'GSE261448', 'Innate immune cells after s.c. BCG vs PBS (PBS controls)', 'NK, Ly6C+ monocytes, neutrophils',
     'Spleen (NK), bone marrow (monocytes, neutrophils)', '3', 'Paired, 2x150', '22-33M pairs',
     'Yes: M1 NK cells and monocytes (PBS controls; neutrophils are not one of the five cell types)', '', 'No',
     'Sex and age not given; each cell type from one tissue.'),
    (M, 'GSE307667', 'Arid2 conditional knockout in B cells (Mb1-Cre controls)', 'Follicular B (also bone-marrow pro-B '
     'and pre-B)', 'Spleen', '4', 'Paired, 2x151', '42-43M pairs', 'Yes: M2 B cells (Mb1-Cre controls)', '', 'No',
     NOGROUPS),
    (M, 'GSE184496', 'Sorted naive CD4 T cells of young and old mice', 'Naive CD4 T', 'Lymph nodes and spleen',
     '6 per age group (3 male, 3 female)', 'Single, 1x100', '12-19M', 'Yes: M2 CD4 T cells (young mice only)', '',
     'Yes: mouse age (7-11 weeks vs 83-86 weeks)',
     'Not used for sex: males and females differ in age (7-8 vs 10-11 weeks).'),
    (M, 'GSE263132', 'Naive CD8 T cells of wild-type C57BL/6J and MD4 mice', 'Naive CD8 T', 'Spleen', '5',
     'Paired, 2x151', '44-49M pairs', 'Yes: M2 CD8 T cells (wild type)', '', 'No',
     'All mice female, 10 weeks; one tissue.'),
    (M, 'GSE236819', 'Ire1/Xbp1 knockout in NK cells (WT littermate controls)', 'NK', 'Spleen', '6',
     'Single, 1x76 (read 1; the SRA runs also hold an index read)', '24-53M', 'Yes: M2 NK cells (WT controls)', '',
     'No', NOGROUPS),
    (M, 'GSE254202', 'C/EBPa enhancer knockout monocytes (wild-type controls)', 'Ly6Chi and Ly6Clo monocytes',
     'Spleen', '3', 'Paired, 2x151', '34-36M pairs', 'Yes: M2 monocytes (wild-type Ly6Chi monocytes)', '', 'No',
     NOGROUPS),
    (M, 'GSE231772', 'Follicular B cells of control and T follicular regulatory cell-deficient mice, naive and '
     'Brucella infected', 'Follicular B', 'Not given', '4 (naive BCL6fl/fl)', 'Paired, 2x100', '29-54M pairs',
     'Yes: M3 B cells (uninfected BCL6fl/fl controls)', '', 'No', 'Sex, age and tissue not given.'),
    (M, 'GSE134322', 'Naive CD4 T cells of wild-type and PIK3CD gain-of-function mice', 'Naive CD4 T', 'Spleen', '3',
     'Paired, 2x150', '23-24M pairs', 'Yes: M3 CD4 T cells (wild-type littermates)', '', 'No',
     'All mice male, 8 weeks; one tissue.'),
    (M, 'GSE284441', 'Aiolos and virtual memory CD8 T cells (naive wild-type CD8 T cells)', 'Naive CD8 T (also '
     'virtual memory CD8 T)', 'Spleen', '4', 'Paired, 2x150', '69-89M pairs', 'Yes: M3 CD8 T cells (naive, wild type)',
     '', 'No', NOGROUPS),
    (M, 'GSE180053', 'METTL3-deficient NK cells (Mettl3fl/fl controls)', 'NK', 'Spleen', '3', 'Paired, 2x150',
     '37-45M pairs', 'Yes: M3 NK cells (Mettl3fl/fl controls)', '', 'No', NOGROUPS),
    (M, 'PRJNA797236', 'Dendritic cell type 3 arises from Ly6C+ monocyte-dendritic cell progenitors (bulk RNA-seq)',
     'Ly6Chi, Ly6Cint and Ly6Clo monocytes, dendritic cells and their progenitors', 'Blood, spleen, bone marrow',
     '4 per tissue (Ly6Chi monocytes, blood and spleen)', 'Paired, 2x150', '15-23M pairs',
     'Yes: M3 monocytes (Ly6Chi monocytes, blood samples only)', '',
     'Yes: mouse tissue (Ly6Chi monocytes, blood vs spleen, 4 vs 4 male mice)', ''),
    # not used in the main analysis
    (M, 'GSE116177', 'Choi et al. 2019; Haemopedia mouse', 'B, CD4 T, CD8 T, NK, monocytes',
     'Spleen and lymph nodes, bone marrow, blood, peritoneal cavity', '2 per cell population',
     'Single 1x100 or paired 2x100', '20-38M', 'No', 'Fewer than three replicates per cell type.',
     'Yes: mouse tissue (Ly6Chi monocytes, bone marrow vs blood)', ''),
    (M, 'ENCODE4: ENCSR968ZOA, ENCSR835KFT', 'Wold lab; B6CASTF1/J F1 hybrid spleen immune cells, total RNA-seq',
     'Naive CD4 T', 'Spleen', '1', 'Single, 1x101', '37-40M', 'No', 'Fewer than three replicates per cell type.', 'No',
     'One female and one male mouse: no replicates per sex.'),
    (M, 'ENCODE4: ENCSR703ACF, ENCSR788PEE, ENCSR243UCH, ENCSR092BBW', 'Wold lab; B6CASTF1/J F1 hybrid',
     'CD19+ B (total B)', 'Spleen', '4', 'Single, 1x101', '~34-44M', 'No',
     'Total, not naive or follicular, B cells; F1 hybrid mice.', 'Yes: mouse sex (8 female vs 8 male mice)', ''),
    (M, 'GSE233347', 'Ly49H+ NK cells of uninfected mice', 'Ly49H+ NK', 'Spleen, liver', '-', 'Paired, 2x101', '-',
     'No', 'Ly49H+ NK cells only (an NK cell subset), not total NK cells.',
     'Yes: mouse tissue (NK cells, spleen vs liver of the same mice)', ''),
    (M, 'GSE112304', 'Quinn et al. 2018; naive CD8 T cells of young and old mice', 'Naive CD8 T', '-',
     '2 pooled samples per age group', 'Paired, 2x100', '-', 'No', 'Fewer than three replicates per cell type.',
     'No', 'Only 2 vs 2 pooled samples.'),
    (M, 'ENCODE4: ENCSR769LJY', 'Wold lab; B6CASTF1/J F1 hybrid', 'Naive CD8 T', 'Spleen', '1', 'Single, 1x101',
     '~34-44M', 'No', 'One replicate only.', 'No', SAME),
    (M, 'ENCODE4: ENCSR959ZHO, ENCSR536HHA', 'Wold lab; B6CASTF1/J F1 hybrid', 'Monocytes', 'Spleen', '1',
     'Single, 1x101', '~34-44M', 'No', 'One replicate per experiment (one male, one female); F1 hybrid mice.',
     'No', 'One replicate per sex.'),
    (M, 'ENCODE4: CD3 T (ENCSR034MGF, ENCSR377UPT), neutrophils (ENCSR251XVY, ENCSR293FRO, ENCSR309GVQ, '
     'ENCSR167IUG, ENCSR698PIS)', 'Wold lab; B6CASTF1/J F1 hybrid', 'Total CD3 T, neutrophils', 'Spleen', '-',
     'Single, 1x101', '~34-44M', 'No', 'Not one of the five cell types (CD3 T cells are not separated into CD4 '
     'and CD8).', 'No', SAME),
    (M, 'GSE180020', 'ImmGen cytokine injections (PBS controls)', 'Follicular B, NK, Ly6C+ monocytes, naive CD4 T, '
     'effector/memory CD8 T', 'Spleen', '9-13', 'Paired, 2x38', '25-43M', 'No',
     'Short reads (2x38) leave many junctions uncovered; the CD8 T cells are effector/memory, not naive.',
     'No', 'Short reads (2x38).'),
    (M, 'GSE122597', 'ImmGen ultra-low-input, deep sequencing', 'Follicular B, NK, naive CD4 T, gamma-delta T, '
     'peritoneal macrophages', 'Spleen, peritoneum', 'One pooled sample per cell type sequenced many times '
     '(technical replicates)', 'Paired, 2x38', '6-10M pairs per library', 'No',
     'No biological replicates; no CD8 T cells or monocytes; short reads (2x38).', 'No',
     'No biological replicates; sex, age and tissue do not vary between the libraries of a cell type.'),
    (M, 'GSE109125', 'ImmGen ultra-low-input, systemwide', '127 immune populations incl. follicular B, naive CD4 T, '
     'naive CD8 T, NK, Ly6C+ monocytes', 'Spleen, blood', '2-4', 'Paired, 2x25 or 2x38', '4-15M pairs', 'No',
     'Very short reads (2x25 or 2x38) and shallow.', 'No', SAME),
    (M, 'GSE124829', 'ImmGen ultra-low-input, sex and age', 'B, CD4 T, CD8 T, NK', 'Spleen', '3 per sex (6 weeks); '
     '1 per sex per age group', 'Paired, 2x25', '4-11M pairs', 'No', 'Short reads (2x25).', 'No',
     'Very short reads (2x25); one mouse per sex per age group.'),
    (M, 'GSE127267, GSE239946, GSE302070, GSE128962', 'ImmGen ultra-low-input (cytokines, IL-2, Klf12)',
     'B, CD4 T, CD8 T, NK, Ly6C+ monocytes (control arms)', 'Spleen', '3-5', 'Paired, 2x25-2x38', '4-11M pairs',
     'No', 'Short reads (2x25-2x38).', 'No', SAME),
    (M, 'GSE116348', 'miR-155 in four immune cell types (wild-type samples)', 'B, CD4 T, dendritic cells, '
     'macrophages', 'Spleen, lymph nodes', '3', 'Paired, 2x50', '22-30M', 'No',
     'All cells activated in vitro (LPS + IL-4, CD3/CD28); short reads; no CD8 T, NK or monocytes.', 'No', SAME),
    (M, 'GSE184841', 'Multiple immune cell types (ILC2 study)', 'CD4 T, CD8 T, B, NK', 'Not given', '2',
     'Single, 1x51', '~15M', 'No', 'Single-end 51 nt reads, shallow; no monocytes.', 'No', SAME),
    (M, 'GSE60101', 'Lara-Astiaso et al. 2014; 16 hematopoietic cell types', 'B, CD4 T, CD8 T, NK, monocytes',
     'Bone marrow, spleen', '1', "Single, ~1x50, 3' tag", '~3M', 'No',
     "3'-end tag sequencing (no splice junctions along the transcript); one replicate.", 'No', SAME),
    (M, 'GSE319684', "Bulk 3' RNA-seq (QuantSeq) of sorted lung cells of SARS-CoV-2 or influenza infected mice",
     'Ly6C+ monocytes, alveolar macrophages, epithelial cells', 'Lung', '4 per group', "3' QuantSeq", '-', 'No',
     "3'-end tag sequencing (no splice junctions along the transcript); lung cells of infected mice.", 'No', SAME),
    (M, 'GSE293874', 'Bulk RNA-seq of fibroblasts reprogrammed into dendritic cells by transcription factors',
     'Transduced fibroblasts, dendritic cells', '-', '-', '-', '-', 'No', 'Not one of the five cell types.',
     'No', SAME),
    (M, 'PRJNA755992', 'Monocyte exhaustion mediated by TRAM (bulk RNA-seq)', 'Monocytes', '-', '1 (PBS)',
     'Paired, 2x75', '~78M pairs', 'No', 'One control replicate; the other samples are LPS treated.', 'No', SAME),
    (M, 'PRJEB55370', 'Bulk RNA-seq of splenic NK cells of Hif1a fl/fl (control) and Ncr1-Cre Hif1a fl/fl mice',
     'NK', 'Spleen', '3 (Hif1a fl/fl)', 'Paired, 2x75', '42-49M pairs', 'No',
     'Three mouse NK cell datasets were sufficient.',
     'No', 'One female and two male control mice (no female replicates); one age group; one tissue.'),
    (M, 'PRJNA531485', 'Avidity selection of the NK response to MCMV (bulk RNA-seq)', 'NK (Ly49H high and low '
     'subsets)', '-', '-', 'Paired, 2x51', '38-45M pairs', 'No', 'NK cell subsets rather than total NK cells.', 'No', SAME),
    (M, 'PRJNA1219535', 'Asparagine deprivation and CD8 T cell activation (bulk RNA-seq)', 'CD8 T (activated)', '-',
     '-', 'Paired, 2x150', '23-27M pairs', 'No', 'Activated cells.', 'No', SAME),
    (M, 'GSE182858', 'Co-housed (CoH) mice, baseline', 'PBMC', 'Peripheral blood', '13', 'Single, 1x151', '-',
     'No', 'Only 17-51% of the reads aligned to the genome (mostly rRNA); PBMC are not a sorted cell type; '
     'the co-housed mice were exposed to pet-store mice.', 'No', SAME),
]

SITE = {'DICE': 'https://dice-database.org', 'BLUEPRINT': 'https://www.blueprint-epigenome.eu',
        'ImmuNexUT': 'https://www.immunexut.org'}


def links(acc):
    """Database links of every accession named in the Accession cell."""
    out = []
    for a in re.findall(r'GSE\d+|E-MTAB-\d+|ENCSR\w+|PRJ[A-Z]+\d+|phs\d+|DICE|BLUEPRINT|ImmuNexUT', acc):
        if a.startswith('E-MTAB'):
            out.append(f'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{a}')
        elif a.startswith('GSE'):
            out.append(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}')
        elif a.startswith('ENCSR'):
            out.append(f'https://www.encodeproject.org/experiments/{a}/')
        elif a.startswith('PRJ'):
            out.append(f'https://www.ncbi.nlm.nih.gov/bioproject/{a}')
        elif a.startswith('phs'):
            out.append(f'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id={a}')
        else:
            out.append(SITE[a])
    return '; '.join(dict.fromkeys(out))


NAME = {'BLUEPRINT': 'BLUEPRINT (EGA)', 'ImmuNexUT': 'ImmuNexUT (NBDC)'}
# biological replicates of each ENCODE experiment (ENCODE portal, November 2026)
ENCODE_REPS = {'ENCSR703ACF': '4 (female)', 'ENCSR788PEE': '4 (female)', 'ENCSR243UCH': '4 (male)',
               'ENCSR092BBW': '4 (male)', 'ENCSR968ZOA': '1 (female)', 'ENCSR835KFT': '1 (male)',
               'ENCSR769LJY': '1 (male)', 'ENCSR959ZHO': '1 (male)', 'ENCSR536HHA': '1 (female)',
               'ENCSR034MGF': '3 (female)', 'ENCSR377UPT': '4 (male)', 'ENCSR251XVY': '2 (male)',
               'ENCSR293FRO': '1 (female)', 'ENCSR309GVQ': '1 (female)', 'ENCSR167IUG': '1 (male)',
               'ENCSR698PIS': '2 (male)'}
REP_COL = 'Biological replicates per cell type'


def add_link_column(df):
    """One row per accession (rows naming several accessions are split), with its database link."""
    rows = []
    for _, r in df.iterrows():
        acc = r['Accession']
        if 'phs' in acc:  # DICE: the raw data are in dbGaP
            toks = re.findall(r'phs\d+', acc)
        else:
            toks = re.findall(r'GSE\d+|E-MTAB-\d+|ENCSR\w+|PRJ[A-Z]+\d+|BLUEPRINT|ImmuNexUT', acc)
        for t in toks:
            new = r.copy()
            if 'phs' in acc:
                pass
            elif t.startswith('ENCSR'):
                new['Accession'] = f'ENCODE4: {t}'
                if REP_COL in new:
                    new[REP_COL] = ENCODE_REPS[t]
                    if new['Cell types'].startswith('Total CD3 T, neutrophils'):
                        new['Cell types'] = 'Total CD3 T' if t in ('ENCSR034MGF', 'ENCSR377UPT') \
                            else 'Neutrophils'
            elif len(toks) > 1 or t in NAME:
                new['Accession'] = NAME.get(t, t)
            new['Link'] = links(t if 'phs' not in acc else t)
            rows.append(new)
    out = pd.DataFrame(rows)
    cols = list(df.columns)
    cols.insert(cols.index('Accession') + 1, 'Link')
    return out[cols].reset_index(drop=True)


def hyperlink(path):
    """Make the cells of the Link column with a single URL clickable."""
    wb = openpyxl.load_workbook(path)
    for ws in wb.worksheets:
        col = [c.column for c in ws[3] if c.value == 'Link'][0]
        for row in ws.iter_rows(min_row=4, min_col=col, max_col=col):
            c = row[0]
            if c.value and ';' not in c.value:
                c.hyperlink = c.value
                c.font = Font(name='Times New Roman', color='0563C1', underline='single')
    wb.save(path)


def main():
    df = add_link_column(pd.DataFrame(ROWS, columns=COLS))
    note = ('Datasets of this study and candidate datasets from the searches for additional public RNA-seq datasets '
            'in GEO, SRA, BioProject, ENA and the ENCODE portal, one row per accession. Main analysis: the GEO and BioProject '
            'queries of Methods 2.2 (human datasets that profile several immune cell types; mouse datasets of the '
            'individual cell types; datasets published until 19 October 2026). Control analysis: search of '
            '19 October 2026. Main analysis: datasets of '
            'sorted or purified immune cells of the five cell types (B, CD4 T, CD8 T, NK cells, monocytes) from unstimulated '
            'cells of healthy donors or unmanipulated or control mice, with at least three biological replicates per '
            'cell type, reads of at least 75 nucleotides and public raw reads. The human datasets (H1, H2, H3) profile '
            'all five cell types; no mouse dataset met these criteria for all five cell types, so the mouse datasets '
            'were additionally searched per cell type (19 October 2026; purified cells of one cell type, at least three replicates of '
            'untreated or control mice), and M1, M2 and M3 were each assembled from one dataset per cell type. '
            'Control analysis: two groups of samples of '
            'one dataset, with replicates, that differ in only sex, age or tissue (Table S13). Reads: per mate for paired-end data.')
    bst.write('Table_S3', 'Supplementary Table S3: Datasets of this study and candidate datasets', note,
              {'datasets': df})
    hyperlink(bst.OUT + 'Table_S3.xlsx')


if __name__ == '__main__':
    main()
