#!/usr/bin/env python3
"""
Supplementary tables of the manuscript (short versions): only the data that the manuscript reports or cites.
The full tables written by build_supp_tables_116.py, confounder_outputs_116.py, table_dataset_search_116.py and
build_sample_metadata_116.py are first moved to supplementary/full_tables/ (done here), and the short tables are
written as supplementary/Table_S<n>.xlsx. Tables S4, S5 and S10 are already short and are kept as they are.
Run after the full tables are (re)built.
"""
import os
import shutil
import sys

import openpyxl
import pandas as pd

import build_supp_tables_116 as bst
import table_numbers_116 as tn
import table_dataset_search_116 as tds

SUP = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_6ds/supplementary/'
FULL = SUP + 'full_tables/'
SHORT = ['S1', 'S2', 'S3', 'S6', 'S8', 'S9', 'S11', 'S12', 'S13']   # internal names (published: table_numbers_116)
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']


def read(n, sheet=0, header=2):
    return pd.read_excel(FULL + f'Table_{tn.name(n)}.xlsx', sheet_name=sheet, header=header)


def main(only=None):
    """only: table names to (re)build, e.g. ['S2', 'S3']; default all."""
    os.makedirs(FULL, exist_ok=True)
    for n in SHORT:   # a freshly built full table replaces the copy in full_tables/
        n = tn.name(n)
        src = SUP + f'Table_{n}.xlsx'
        if os.path.exists(src) and (not os.path.exists(FULL + f'Table_{n}.xlsx') or
                                    os.path.getmtime(src) > os.path.getmtime(FULL + f'Table_{n}.xlsx') + 1):
            if openpyxl.load_workbook(src, read_only=True).properties.keywords != 'short table':  # a full table
                shutil.copy2(src, FULL + f'Table_{n}.xlsx')
    bst.OUT = SUP
    write = bst.write

    def tagged(name, *a, **k):   # mark the short tables, so that a rerun does not copy them to full_tables/
        write(name, *a, **k)
        name = tn.name(name)
        wb = openpyxl.load_workbook(SUP + f'{name}.xlsx')
        wb.properties.keywords = 'short table'
        wb.save(SUP + f'{name}.xlsx')
    bst.write = tagged

    def s1():
        # S1: orthology counts per definition and exon-rank consistency (Results 3.1)
        s1 = read('S1', 'orthology_map')[['variant', 'definition', 'genes_with_ortholog_exons', 'exons_h',
                                          'exons_h_with_ortholog', 'exons_m', 'exons_m_with_ortholog', 'unique_boundaries']]
        rc = read('S1', 'rank_consistency')[['orthology', 'in_canonical_transcript_of_both', 'pct_same_rank',
                                            'pct_constant_offset_within_gene']]
        bst.write('Table_S1', 'Supplementary Table S1: Human-mouse exon orthology',
                  {'orthology': 'For each definition of exon orthology: genes with orthologous exons, human (h) and mouse '
                                '(m) exons and those with an ortholog, and unique orthologous exon boundaries.',
                   'rank_consistency': 'For the exon boundaries used by the alternative splicing events, in canonical '
                                       'transcripts of both species: % with the same exon rank in human and mouse, and % '
                                       'with a constant rank offset within the gene.'},
                  {'orthology': s1, 'rank_consistency': rc})

    def s2():
        # S2: one row per sample, the fields named in Methods 2.2
        s2 = read('S2', header=1)
        keep = ['dataset', 'species', 'GSM', 'sample_name_used', 'group_name', 'tissue_source',
                'sorting_definition_markers', 'cell_state', 'treatment_condition', 'donor_or_mouse_id', 'sex', 'age',
                'strain_or_ethnicity', 'read_length_nt', 'flagstat_mapped']
        bst.write('Table_S2', 'Supplementary Table S2: Samples',
                  'Samples used: dataset, accession, cell type group, tissue, sorting markers, cell state, treatment, donor or '
                  'mouse identity, sex, age, strain, read length and mapped reads.', {'samples': s2[[c for c in keep if c in s2.columns]]})

    def s3():
        # S3: dataset, whether used and why not
        s3 = read('S3')[['Species', 'Accession', 'Link', 'Cell types', 'Main analysis (human vs mouse)',
                         'Reason not used in the main analysis', 'Control analysis (sex, age, tissue)',
                         'Reason not used in the control analysis']]
        bst.write('Table_S3', 'Supplementary Table S3: Datasets of this study and candidate datasets',
                  'Main analysis: the GEO and BioProject queries of Methods 2.2 (published until 19 October 2026); control '
                  'analysis: search of 19 October 2026. Selection criteria in Methods 2.2.',
                  {'datasets': s3})
        tds.hyperlink(SUP + 'Table_S3.xlsx')   # clickable links, as in the full table

    def s6():
        # S6: only the sensitivity results reported in the text
        orth = read('S6', 'orthology_cutoff')[['orthology', 'cell', 'testable', 'pct_DS', 'DS events (any cell type)']]
        ce = read('S6', 'common_events')
        ss = read('S6', 'splice_site_subset')
        ss = ss[ss['subset'].isin(['all', 'correspondence'])][['subset', 'cell', 'testable', 'pct_DS']]
        mu = read('S6', 'minor_isoform_usage')[['usage_bin', 'DS', 'splicing conserved', 'inconclusive', 'DS %', 'splicing conserved %']]
        sh = read('S6', 'small_shift_A3SS_A5SS')[['event_type', 'pct_DS_small', 'pct_DS_larger', 'fisher_p']]
        sc = read('S6', 'single_cell_type_DS')[['DS_in_one_cell_type', 'status_in_other_cell_types', 'percent']]
        sp = read('S6', 'species_vs_mixed_splits')
        bst.write('Table_S6', 'Supplementary Table S6: Sensitivity analyses', {
            'orthology_cutoff': 'Human-mouse classification with each definition of exon orthology (id50, id70, id80: at '
                                'least 50%, 70% or 80% identity; id50_rbh: 50% identity and reciprocal best match): '
                                'testable events and % differentially spliced (DS) per immune cell type, and the number '
                                'of events DS in at least one cell type.',
            'common_events': 'Events testable in all 15 dataset pairs of a cell type (common_events): for each pair, '
                             'the median max|dPSI| of the events and the Pearson r of the junction PSIs (Figure 3). '
                             'kind: HH, human-human; MM, mouse-mouse; cross, human-mouse.',
            'splice_site': 'Testable events and % differentially spliced per cell type, among all events and among the '
                           'events with well-supported human-mouse splice-site correspondence (Methods).',
            'minor_isoform': 'Human-mouse classifications (event x cell type) by the usage of the minor isoform (PSI of '
                             'the less used isoform, averaged over the six datasets): number and % differentially '
                             'spliced, splicing conserved and inconclusive.',
            'small_shift': "Alternative 3' (A3SS) and 5' (A5SS) splice site events: % differentially spliced among "
                           'events with alternative splice sites up to 6 nt apart (small) and more distant (larger), '
                           "and Fisher's exact test P.",
            'single_cell_type_DS': 'Events differentially spliced in only one immune cell type: their status in the '
                                   'other cell types (% of the event x cell type combinations).',
            'species_vs_mixed_splits': 'Number of differentially spliced events per cell type when the classification '
                                       'is applied to each split of the six datasets into two groups of three: the '
                                       'human-mouse split and the nine splits that mix the species.'},
                  {'orthology_cutoff': orth, 'common_events': ce, 'splice_site': ss,
                   'minor_isoform': mu, 'small_shift': sh, 'single_cell_type_DS': sc,
                   'species_vs_mixed_splits': sp})

    def s7():
        # S7: odds ratios and P values
        s7 = read('S7')[['outcome', 'term', 'odds_ratio', 'p']]
        bst.write('Table_S7', 'Supplementary Table S7: Effect of coverage on testability and classification',
                  'Logistic regression over event x cell type (cell type as covariate): odds ratio and P value of event '
                  'coverage (log10 reads per sample), number of junctions and number of events per gene.', {'model': s7})

    def s8():
        # S8: the orthologous coordinates of every junction (cited for the mouse coordinates)
        x = pd.ExcelFile(FULL + f'Table_{tn.name("S8")}.xlsx')
        j = pd.concat([pd.read_excel(x, s, header=2, usecols=range(8)) for s in x.sheet_names]).drop_duplicates(
            'Human junction')
        bst.write('Table_S8', 'Supplementary Table S8: Orthologous junctions',
                  'Each junction of the alternative splicing events with its orthologous mouse junction, exon ranks, event, '
                  'gene and event type.', {'junctions': j})

    def s9():
        # S9: the single-pair results shown in Figure 2A-C
        s9 = read('S9')[['comparison', 'kind', 'cell', 'testable events', '% differentially spliced',
                         '% splicing conserved']]
        s9['% inconclusive'] = (100 - s9['% differentially spliced'] - s9['% splicing conserved']).round(1)
        bst.write('Table_S9', 'Supplementary Table S9: Single dataset pairs',
                  'Per dataset pair and cell type: testable events and % differentially spliced (LeafCutter adjusted '
                  'P <= 0.05 and max|dPSI| >= 0.1), splicing conserved (LeafCutter P > 0.1 and max|dPSI| <= 0.05) and inconclusive '
                  '(Figure 2A-C).',
                  {'comparisons': s9})

    def s11():
        # S11: the class of each event in each cell type and in fibroblasts
        s11 = read('S11')[['Alternative splicing event', 'Gene', 'Type'] + CELLS + ['Fibroblasts (human vs mouse)']]
        bst.write('Table_S11', 'Supplementary Table S11: Classification of each event',
                  'Per immune cell type: differentially spliced (DS), splicing conserved, inconclusive or not testable (Table S10, '
                  'Methods 2.5); fibroblasts: single human-mouse pair.', {'events': s11})

    def s12():
        # S12: the significant terms
        x = pd.ExcelFile(FULL + f'Table_{tn.name("S12")}.xlsx')
        sig = []
        for s in x.sheet_names:
            if s == 'coverage_check':
                continue
            d = pd.read_excel(x, s, header=2)
            d = d[d['Adjusted P-value'] <= 0.1][['Term', 'Overlap', 'Adjusted P-value']]
            sig.append(d.assign(list=s.split('_')[0], library=s.split('_', 2)[-1]))
        sig = pd.concat(sig)[['list', 'library', 'Term', 'Overlap', 'Adjusted P-value']]
        bst.write('Table_S12', 'Supplementary Table S12: Functional enrichment',
                  'Significant terms (FDR <= 0.1) of the genes with only splicing conserved events and of the genes with '
                  'only differentially spliced events.',
                  {'significant_terms': sig})

    def s13():
        # S13: the control results reported in Supplementary Note 1
        s13 = read('S13')[['comparison', 'dataset', 'cell_type', 'common_events', 'pct_DS_control',
                           'median_abs_dpsi_control', 'median_abs_dpsi_cross', 'pearson_control', 'pearson_cross']]
        bst.write('Table_S13', 'Supplementary Table S13: Sex, age and tissue controls',
                  'Per control comparison: events (testable in the control and in all nine human-mouse pairs), % '
                  'differentially spliced in the control, median max|dPSI| and Pearson r in the control and in the '
                  'human-mouse pairs (cross; mean over the nine pairs).', {'controls': s13})

    steps = {'S1': s1, 'S2': s2, 'S3': s3, 'S6': s6, 'S8': s8, 'S9': s9, 'S11': s11, 'S12': s12, 'S13': s13}
    for n, step in steps.items():
        if only and n not in only:
            continue
        if not os.path.exists(FULL + f'Table_{tn.name(n)}.xlsx'):   # full table not built yet (e.g. a partial rerun)
            print(f'skipped Table_{n}: no full table in', FULL)
            continue
        step()


if __name__ == '__main__':
    main(sys.argv[1:] or None)
