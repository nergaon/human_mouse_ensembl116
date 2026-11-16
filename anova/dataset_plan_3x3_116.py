#!/usr/bin/env python3
"""Datasets that can give 3 human and 3 mouse datasets per immune cell type, each with >= 3 biological replicates.
Facts from Table S2 (datasets in use), ensembl116/m3/ (M3 alignments) and ensembl116/metadata/m3_candidates.tsv.
Output: ensembl116/m3/dataset_plan_3x3.xlsx"""
import pandas as pd
import build_supp_tables_116 as bst

H, M = 'Human', 'Mouse'
IN, DONE, TODO = 'in the study', 'downloaded and aligned', 'to download'
R = [
 # cell, species, dataset, replicates, reads, status, note, proposed
 ('B', H, 'GSE60424 (H1)', 4, 'paired 2x50', IN, '', 'yes'),
 ('B', H, 'GSE115736 (H2)', 5, 'single 1x100', IN, '', 'yes'),
 ('B', H, 'GSE118165 (H3)', 4, 'paired ~2x98', IN, 'naive B', 'yes'),
 ('B', H, 'GSE107011', 4, 'paired 2x50', DONE, 'naive B; alternative (sex control dataset)', 'alternative'),
 ('B', M, 'GSE287527 (M2)', 4, 'paired 2x101', IN, 'follicular B', 'yes'),
 ('B', M, 'GSE307667', 4, 'paired 2x151', DONE, 'follicular B; Mb1-Cre controls (Cd79a+/-)', 'yes'),
 ('B', M, 'GSE231772', 4, 'paired 2x100', DONE, 'follicular B; Bcl6 fl/fl controls', 'yes'),
 ('B', M, 'ENCODE4 B (4 experiments)', 16, 'single 1x101', DONE, 'total CD19+ B, F1 hybrid; sex control', 'alternative'),
 ('B', M, 'GSE116177 (M1)', 2, 'single 1x100', IN, 'only 2 replicates', 'no'),
 ('CD4 T', H, 'GSE60424 (H1)', 4, 'paired 2x50', IN, 'total CD4 T', 'yes'),
 ('CD4 T', H, 'GSE115736 (H2)', 5, 'single 1x100', IN, 'total CD4 T', 'yes'),
 ('CD4 T', H, 'GSE118165 (H3)', 4, 'paired ~2x98', IN, 'naive CD4 T', 'yes'),
 ('CD4 T', M, 'GSE215775', 3, 'paired 2x151', DONE, 'naive CD4 T, spleen', 'yes'),
 ('CD4 T', M, 'GSE184496 (young)', 6, 'single 1x100', DONE, 'naive CD4 T, 7-11 weeks; age control dataset', 'yes'),
 ('CD4 T', M, 'GSE134322', 3, 'paired 2x150', TODO, 'naive CD4 T, negatively isolated, male 8 weeks', 'yes'),
 ('CD4 T', M, 'GSE116177 (M1)', 8, 'single 1x100', IN, '2 per subset (naive, total, effector, memory)', 'no'),
 ('CD4 T', M, 'ENCODE4 CD4 (M2)', 2, 'single 1x101', IN, 'only 2 replicates', 'no'),
 ('CD8 T', H, 'GSE60424 (H1)', 4, 'paired 2x50', IN, 'total CD8 T', 'yes'),
 ('CD8 T', H, 'GSE115736 (H2)', 5, 'single 1x100', IN, 'total CD8 T', 'yes'),
 ('CD8 T', H, 'GSE118165 (H3)', 4, 'paired ~2x98', IN, 'naive CD8 T', 'yes'),
 ('CD8 T', M, 'GSE313371 (M2)', 3, 'paired 2x150', IN, 'CD8 T, uninfected wild type', 'yes'),
 ('CD8 T', M, 'GSE263132', 5, 'paired 2x151', DONE, 'naive CD8 T, CD44-low', 'yes'),
 ('CD8 T', M, 'GSE284441', 4, 'paired 2x150', TODO, 'naive CD8 T (bead-isolated)', 'yes'),
 ('CD8 T', M, 'GSE333264 / GSE305926', '3 / 3+3', 'paired 2x150', TODO, 'naive CD8 T; alternatives', 'alternative'),
 ('CD8 T', M, 'GSE116177 (M1)', 6, 'single 1x100', IN, '2 per subset (naive, total, memory)', 'no'),
 ('NK', H, 'GSE60424 (H1)', 4, 'paired 2x50', IN, '', 'yes'),
 ('NK', H, 'GSE115736 (H2)', 5, 'single 1x100', IN, '', 'yes'),
 ('NK', H, 'GSE118165 (H3)', 4, 'paired ~2x98', IN, 'mature NK', 'yes'),
 ('NK', M, 'GSE261448 (M2)', 3, 'paired 2x150', IN, 'spleen NK, PBS controls', 'yes'),
 ('NK', M, 'GSE236819', 6, 'single 76', DONE, 'spleen NK, steady state', 'yes'),
 ('NK', M, 'GSE180053', 3, 'paired 2x150', DONE, 'spleen NK, Mettl3 fl/fl controls', 'yes'),
 ('NK', M, 'GSE305179 / GSE342137', '3 / 4', 'paired 2x100 / 2x150', TODO, 'alternatives', 'alternative'),
 ('NK', M, 'GSE116177 (M1)', 2, 'paired 2x100', IN, 'only 2 replicates', 'no'),
 ('Monocytes', H, 'GSE60424 (H1)', 4, 'paired 2x50', IN, '', 'yes'),
 ('Monocytes', H, 'GSE115736 (H2)', 5, 'single 1x100', IN, '', 'yes'),
 ('Monocytes', H, 'GSE118165 (H3)', 4, 'paired ~2x98', IN, '', 'yes'),
 ('Monocytes', M, 'GSE261448 (M2)', 3, 'paired 2x150', IN, 'Ly6C+ bone marrow, PBS controls', 'yes'),
 ('Monocytes', M, 'GSE254202', 3, 'paired 2x151', DONE, 'Ly6Chi spleen', 'yes'),
 ('Monocytes', M, 'PRJNA797236 (blood)', 4, 'paired 2x150', DONE, 'Ly6Chi blood; tissue control dataset', 'yes'),
 ('Monocytes', M, 'GSE116177 (M1)', 6, 'single 1x100', IN, '2 per population (blood, bone marrow, inflammatory)', 'no'),
]
t = pd.DataFrame(R, columns=['Cell type', 'Species', 'Dataset', 'Biological replicates', 'Reads', 'Status', 'Notes',
                             'Proposed for 3 x 3'])
bst.OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/m3/'
bst.write('dataset_plan_3x3', 'Datasets for 3 human and 3 mouse datasets per cell type (>= 3 replicates each)',
          'Proposed = three human and three mouse datasets per cell type, each with at least three biological '
          'replicates of that cell type; alternative = usable instead; no = fewer than three replicates of one '
          'population. Status: in the study, downloaded and aligned (M3 or control datasets), or to download.',
          {'datasets': t})
