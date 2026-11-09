#!/usr/bin/env python3
"""
Supplementary Material document of the revision (replaces the submitted Table_1.DOCX):
Tables S3, S4 and S9 as tables, the list of supplementary data sheets (Excel), and the legends of
Figures S1-S10.  Output: ensembl116/revision_separation/supplementary/Supplementary_Material.docx
"""
import glob
import re

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

from docx_style_116 import apply_style, finalize

B = '/gpfs0/tals/projects/Analysis/human_mouse_exons/'
ED = B + 'ensembl116/events/id50_sep/'
OUT = B + 'ensembl116/revision_separation/supplementary/Supplementary_Material.docx'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_NAME = {'CD4T': 'CD4T', 'CD8T': 'CD8T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mono', 'Fibro': 'Fibroblasts'}
DATASETS = [('H1', 'Human', 'GSE60424'), ('H2', 'Human', 'GSE115736'), ('H3', 'Human', 'GSE118165'),
            ('M1', 'Mouse', 'GSE116177'), ('M2', 'Mouse', 'GSE287527 / GSE261448 / ENCODE4 / GSE313371'),
            ('FH', 'Human', 'E-MTAB-5919'), ('FM', 'Mouse', 'E-MTAB-5919')]

DATA_SHEETS = [
    ('S2', 'Sample-level metadata of all samples used: dataset, GEO/ENCODE and SRA accessions, cell type, sorting '
           'markers, tissue source, cell state, treatment, donor/mouse identity, sex, age, strain, library type, read '
           'length and mapping statistics.'),
    ('S1', 'Human-mouse exon orthology map (Ensembl 116) with orthologous exons defined by at least 50% (main), 70% '
           'or 80% identity, or by reciprocal best matches.'),
    ('S7', 'Junction-level LeafCutter differential splicing results, one sheet per comparison: human and orthologous '
           'mouse junction, human and mouse exon ranks at the junction ends, alternative splicing event, human gene symbol, alternative splicing type, well-supported '
           'splice-site correspondence, and fitted PSI per group, ΔPSI and p.adjust in each cell type.'),
    ('S8', 'Differential splicing output of every comparison and cell type: testable events, % differentially '
           'spliced, % equivalent at margins 0.025, 0.05 and 0.1, median max(|ΔPSI|) and PSI correlation; the '
           'human-mouse pairs (Figure 2C) are a step before the classification of the events from all five datasets '
           '(Table 1).'),
    ('S10', 'Classification of each orthologous alternative splicing event in each immune cell type (confidence '
           'level) and in fibroblasts.'),
    ('S11', 'Functional enrichment of genes with high confidence differentially spliced or equivalent events '
            '(GO biological process 2025, Reactome 2022, hallmark), and coverage-adjusted check of enriched terms.'),
    ('S5', 'Sensitivity analyses: orthology stringency, ΔPSI threshold, events testable in all comparisons, '
            'events with well-supported splice-site correspondence, matched cell populations, classification by '
            'minor-isoform usage, small alternative splice-site shifts, and the status in the other cell types of '
            'events that are differentially spliced with high confidence in only one cell type.'),
    ('S6', 'Effect of coverage, number of junctions and number of events per gene on testability and '
            'classification (logistic regression).'),
    ('S12', 'Sex, age and tissue control comparisons within each species: datasets, compared samples, testable and '
            'differentially spliced events, and comparison with the human-mouse and within-species dataset pairs of '
            'the same cell type on the same events.'),
]

FIG_LEGENDS = [
    ('Figure S1', 'The status of an orthologous alternative splicing event in a single comparison of two datasets, '
                   'or of two sample groups, in one cell type. A. Differential splicing was analyzed in each immune cell '
                   'type between (i) each of the three human and each of the two mouse datasets (used only for effect sizes, '
                   'Figure 3), (ii) each pair of human datasets and (iii) the two mouse datasets, and (iv) between human '
                   'and mouse fibroblasts. B. The comparisons classified in this way (human-mouse events in immune cells '
                   'are classified only from all five datasets, Figure 1B, C), including the sex, age and tissue control comparisons, '
                   'and where their results are shown. C. The decision process: an event is not testable if fewer than the '
                   'minimal number of samples (g) in a group have at least 20 JSRs in the event; differentially '
                   'spliced if the LeafCutter Benjamini-Hochberg adjusted P value (p.adjust) is at most 0.05 and '
                   'max(|ΔPSI|) is at least 0.1; equivalent if the 90% CI of the ΔPSI of every junction lies within '
                   '±0.1; and inconclusive otherwise. The human-mouse classification of each event across all five '
                   'immune cell datasets is shown in Figure 1B, C.'),
    ('Figure S2', 'Classification of human-mouse orthologous alternative splicing events including low confidence '
                  'events. A-C. As in Figure 2, but as percentages of all 4,570 orthologous alternative splicing events, '
                  'including those that are not testable in the comparison (fewer than the minimal number of samples with at '
                  'least 20 JSRs in a group; Figure S1C). D. As in '
                  'Figure 2. E. UpSet plot of differentially spliced events (high or low confidence) '
                  'shared between cell types; as in Figure 2E, divergence in a single cell type does not imply cell type '
                  'specificity (Table S5). F. UpSet plot of equivalent events (high or low confidence). G. Heatmap of '
                  'events that are differentially spliced (green) or equivalent (blue) with high confidence in exactly two '
                  'immune cell types (events with high confidence in three or more are shown in Figure 2G), colored as in '
                  'Figure 2G (dotted, high confidence; plain, low confidence), and their status in fibroblasts (Fb). '
                  'Gene symbols are shown for the containing genes.'),
    ('Figure S3', 'Splicing differences of the classified events in every comparison. For each immune cell type, '
                  'the max(|ΔPSI|) of (A) the differentially spliced and (B) the equivalent events (high or low confidence; '
                  'n, number of events) in each of the three human-human (blue), one mouse-mouse (green) and six '
                  'human-mouse (orange) comparisons in which the event is testable. Dotted line, ΔPSI = 0.1; numbers '
                  'above the boxes, median max(|ΔPSI|). Boxes, interquartile range; line, median; whiskers, 5th-95th '
                  'percentiles.'),
    ('Figure S4', 'The orthologous skipped exon 5 event in PUF60 gene displays differential splicing between human '
                  'and mouse. A. The schematics of the orthologous skipped exon event. Exons are represented by '
                  'rectangles, splicing junctions by arcs. Exon numbers are those of the human Ensembl canonical transcript '
                  'and of the mouse transcript most similar to it that contains the orthologous junctions (Puf60-202, '
                  'ENSMUST00000100527; the mouse Ensembl canonical transcript lacks exons 2-5). B. Stacked bars display the number of junction spanning '
                  'reads that cover the alternative splicing junctions in each sample (top), and the same numbers '
                  'normalized to the total to display PSI values (bottom). Colors are as in A. M2, M1 – mouse '
                  'datasets; H1, H2, H3 – human datasets; M, H – mouse and human fibroblasts.'),
    ('Figure S5', 'The orthologous alternative 3\' splicing event in RPL22L1 exon 2 is differentially spliced '
                  'between human and mouse. A. An orthologous alternative 3\' splicing event of exon 2. The '
                  'transcript that includes the blue junction codes for an additional Glutamine, in the human '
                  'canonical transcript, compared to the transcript that includes the light blue junction, which is '
                  'shorter by three nucleotides. B. As in Figure S4B. Colors are as in A.'),
    ('Figure S8', 'The orthologous alternative first exon in DOCK8 gene displays a similar splicing ratio in human '
                  'and mouse but differential splicing between immune cells. A. An orthologous alternative first '
                  'exon. The transcript that includes the blue junction codes for a longer transcript, compared to '
                  'the transcript with the light blue junction. The mouse transcript with the light blue junction '
                  'biotype is protein coding, CDS not defined. B. As in Figure S4B. Colors are as in A.'),
    ('Figure S9', 'The orthologous complex skipped exons event in PTPRC displays differential splicing between '
                  'immune cell types and species. A. The schematics of the orthologous skipped exons event. Exons are '
                  'represented by rectangles, splicing junctions by arcs. B. As in Figure S4B. Colors are as in A; two '
                  'junctions that are part of the event in the current analysis but are not drawn in A are shown '
                  'in grey (chr1:198,692,373-198,702,387) and purple (chr1:198,699,704-198,704,472).'),
    ('Figure S10', 'The orthologous complex skipped exons event in CD47 displays differential splicing between '
                  'species, but is similar between immune cell types. A. The schematics of the orthologous skipped '
                  'exons event. Exons are represented by rectangles, splicing junctions by arcs. B. As in Figure S4B. '
                  'Colors are as in A; the junction shown in purple in B (chr3:108,055,566-108,057,477) is part of the '
                  'event in the current analysis but is not drawn in A.'),
    ('Figure S6', 'The orthologous skipped exon event in UCKL1 displays differential splicing between human and '
                  'mouse immune cells, and a similar use ratio in human and mouse fibroblasts. A. The schematics of '
                  'the orthologous skipped exon event. Exons are represented by rectangles, splicing junctions by '
                  'arcs. B. As in Figure S4B. Colors are as in A.'),
    ('Figure S7', 'Sex, age and tissue of origin cause much smaller splicing differences than species. Each '
                  'control comparison is between two groups of samples of the same dataset (same cell type, '
                  'laboratory and protocol), within one species (H, human; M, mouse), that differ in one factor: sex (H: GSE107011, male vs female donors; '
                  'M: ENCODE spleen B cells, male vs female mice), age (H: GSE143213, young vs old donors of naive CD4 T '
                  'cells; M: GSE184496, young vs old mice, naive CD4 T cells) and tissue (H: GSE133383, NK cells from '
                  'blood vs spleen of the same donors; M: GSE233347, NK cells from spleen vs liver; GSE116177, '
                  'monocytes from bone marrow vs blood), compared with the human-human, mouse-mouse and human-mouse '
                  'dataset pairs of the same cell type (median over pairs), all tested as single comparisons (Figure S1C). '
                  'A. Percentage of differentially spliced events (the denominator includes the inconclusive and '
                  'equivalent events): for the control comparisons, of the events testable in the control comparison '
                  'and in all six human-mouse pairs of the cell type; for the human-human and mouse-mouse dataset '
                  'pairs, of all events testable in each pair, as in Table 1 (human-mouse events are classified only '
                  'from all five datasets, Table 1). B, C. On the events testable in the control comparison and in all six '
                  'human-mouse pairs of the cell type. B. Median max(|ΔPSI|) of the events. C. Pearson correlation of '
                  'junction PSI. T4, T8 – CD4 and CD8 T cells; B – B cells; NK – NK cells; Mo – monocytes.'),
]


def set_cell_borders(table):
    tbl = table._tbl
    borders = OxmlElement('w:tblBorders')
    for edge in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        el = OxmlElement(f'w:{edge}')
        el.set(qn('w:val'), 'single')
        el.set(qn('w:sz'), '4')
        el.set(qn('w:color'), '000000')
        borders.append(el)
    # schema order in tblPr: ... tblW, jc, tblCellSpacing, tblInd, tblBorders, shd, tblLayout, tblCellMar, tblLook
    later = [c for c in tbl.tblPr if c.tag in {qn(f'w:{t}') for t in ('shd', 'tblLayout', 'tblCellMar', 'tblLook',
                                                                       'tblCaption', 'tblDescription')}]
    if later:
        later[0].addprevious(borders)
    else:
        tbl.tblPr.append(borders)


def add_table(doc, header, rows, font=12, widths=None):
    t = doc.add_table(rows=1, cols=len(header))
    set_cell_borders(t)
    t.autofit = widths is None
    for c, h in zip(t.rows[0].cells, header):
        c.text = ''
        r = c.paragraphs[0].add_run(str(h))
        r.bold = True
        r.font.size = Pt(font)
    for row in rows:
        cells = t.add_row().cells
        for c, v in zip(cells, row):
            c.text = ''
            c.paragraphs[0].add_run(str(v)).font.size = Pt(font)
    for r in t.rows:  # house style: spacing after 0 pt, single line spacing
        for c in r.cells:
            for par in c.paragraphs:
                par.paragraph_format.space_after = Pt(0)
                par.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    if widths:
        for r in t.rows:
            for c, w in zip(r.cells, widths):
                c.width = Cm(w)
    return t


def table_s1():
    design = pd.read_csv(ED + 'leafcutter/sample_design.tsv', sep='\t')
    meta = pd.read_csv(B + 'ensembl116/metadata/sample_metadata_revision.tsv', sep='\t')
    meta['mapped'] = pd.to_numeric(meta['flagstat_mapped'].astype(str).str.replace(',', '').str.extract(r'(\d+)')[0],
                                   errors='coerce')
    # sample columns of the count table are '<sample>_<study>'; E-MTAB-5919 studies are EMTAB5919H / EMTAB5919M
    study_acc = {'EMTAB5919H': 'E-MTAB-5919', 'EMTAB5919M': 'E-MTAB-5919'}
    mapped = {(r['dataset'], str(r['sample_name_used'])): r['mapped'] for _, r in meta.iterrows()}
    design['study'] = design['sample'].str.rsplit('_', n=1).str[1]
    design['name'] = design['sample'].str.rsplit('_', n=1).str[0]
    rows = []
    for label, sp, _ in DATASETS:
        d = design[design['dataset'] == label]
        for cell in CELLS + ['Fibro']:
            for study, s in d[d['cell'] == cell].groupby('study', sort=False):
                acc = study_acc.get(study, study)
                vals = [mapped.get((acc, n)) for n in s['name']]
                vals = [v for v in vals if pd.notna(v)]
                rows.append([f'{label} ({sp})', acc, CELL_NAME[cell], len(s),
                             f'{sum(vals) / len(vals):,.0f}' if vals else 'NA'])
    return rows


def compared_samples():
    """Table S4: one row per comparison (group 1 vs group 2), taken from the LeafCutter comparisons that were run.
    Group names as in the group_name column of Table S2, which lists the samples of each group."""
    MATCHED = {'CD4T', 'CD8T', 'Mono'}  # the matched M1 subset differs from M1 only for these cell types
    order_lab = ['H1', 'H2', 'H3', 'M1', 'M1-matched', 'M2', 'FH', 'FM']
    order_cell = CELLS + ['Fibro']
    pairs = set()
    for line in open(ED + 'leafcutter/jobs.txt'):
        sample_set, pair, cell = line.strip().rstrip('/').split('/')[-3:]
        labs = [lab + ('-matched' if sample_set == 'matched' and lab == 'M1' and cell in MATCHED else '')
                for lab in pair.split('_')]
        pairs.add((cell, *sorted(labs, key=order_lab.index)))
    subset = {'CD4T': 'naive CD4 T cells', 'CD8T': 'naive CD8 T cells', 'Mono': 'blood monocytes'}

    def reason(cell, a, b):
        if cell == 'Fibro':
            return 'Human vs mouse in a non-immune cell type, both species from the same study'
        if a[0] == b[0]:
            sp = 'human' if a[0] == 'H' else 'mouse'
            extra = ' (for the matched analysis)' if 'matched' in a + b else ''
            who = 'donors' if sp == 'human' else 'mice'
            return f'Within-species baseline{extra}: differences between {sp} datasets (laboratory, protocol, {who})'
        if 'matched' in b:
            return (f'Human vs mouse with more similar populations: M1 restricted to {subset[cell]} (M1 also '
                    'includes other cell states or tissues)')
        return 'Human vs mouse (main analysis)'
    return [[f'{a}:{CELL_NAME[cell]}', f'{b}:{CELL_NAME[cell]}', reason(cell, a, b)] for cell, a, b in
            sorted(pairs, key=lambda k: (order_cell.index(k[0]), k[1][0] != k[2][0],
                                         order_lab.index(k[1]), order_lab.index(k[2])))]


def main():
    doc = Document()
    apply_style(doc)
    h = doc.add_heading('Supplementary Material', level=1)
    doc.add_paragraph('Limited conservation of alternative splicing use ratio between human and mouse immune cells. '
                      'Hadas Ner-Gaon, Naama Nudelman, Tal Shay')
    doc.add_heading('Supplementary Tables', level=2)
    doc.add_paragraph('Tables S3, S4 and S9 are given below; all other tables are provided as separate data sheets '
                      '(Excel files named Table_S<number>.xlsx).')
    sheets = dict(DATA_SHEETS)

    def text_s3():
        add_table(doc, ['Dataset', 'Accession', 'Cell type', '#Replicates', 'Average mapped reads'], table_s1())
        doc.add_paragraph('-i (min_samples_per_intron) and -g (min_samples_per_group) of LeafCutter were set to 60% '
                          'of the smallest group, and at least 2: 3 for H2 vs M1 in CD4T, CD8T and Mono, and 2 in all '
                          'other comparisons. Mapped reads: SAMtools flagstat mapped reads (alignments of paired reads '
                          'are counted per mate).')

    def text_s4():
        add_table(doc, ['Group 1', 'Group 2', 'Reason for the comparison'], compared_samples(),
                  widths=[3.5, 3.5, 9])
        doc.add_paragraph('Each row is one comparison of group 1 vs group 2. Group = dataset:cell type (datasets as in '
                          'Table S3); the samples of each group are listed in Table S2 (column group_name). '
                          'M1-matched: the subset of M1 (GSE116177) samples of naive CD4 or CD8 T cells, or of blood '
                          'monocytes (for B and NK cells M1 is used as is). FH vs FM (fibroblasts) were compared '
                          'without confidence levels.')

    def text_s10():
        hc = pd.read_csv(ED + 'classification/summary_by_cell.tsv', sep='\t', index_col=0)
        st = pd.read_csv(ED + 'classification/event_status_by_cell.tsv', sep='\t')
        tc = st[st['status_high'] != 'not testable'].groupby('cell').size()
        rows = [[CELL_NAME[c], tc[c]] + [f'{hc.loc[c, k]} ({100 * hc.loc[c, k] / tc[c]:.1f}%)' for k in
                                         ('high_divergent', 'high_conserved', 'high_low_divergent',
                                          'high_low_conserved')] for c in CELLS]
        add_table(doc, ['Cell type', 'Testable events', 'DS, high', 'Equivalent, high',
                        'DS, high + low', 'Equivalent, high + low'], rows)
        doc.add_paragraph('Testable: every junction of the event has a PSI (mean over the samples with at least 20 '
                          'junction spanning reads in the cluster) in each of the three human and two mouse datasets. '
                          'Differentially spliced (DS): for at least one junction, the PSI of every human dataset differs '
                          'from the PSI of every mouse dataset by at least the gap (high confidence, gap ≥ 0.2; low '
                          'confidence, 0.1 ≤ gap < 0.2). Equivalent: the PSI of every junction lies within a range of at '
                          'most 0.1 over all five datasets (high confidence, at most 0.05; low confidence, more than 0.05). '
                          'High + low confidence includes the high-confidence events. Also provided as Table_S9.xlsx.')

    text_tables = {'S3': ('Datasets summary', text_s3), 'S4': ('Compared samples', text_s4),
                   'S9': ('Differentially spliced and equivalent events by confidence level', text_s10)}
    for n in range(1, 13):
        key = f'S{n}'
        p = doc.add_paragraph()
        if key in text_tables:
            title, fn = text_tables[key]
            p.add_run(f'Table {key}: {title}').bold = True
            fn()
        else:
            p.add_run(f'Table {key} (data sheet, Table_{key}.xlsx): ').bold = True
            p.add_run(sheets[key])

    doc.add_heading('Supplementary Figures Legends', level=2)
    for name, text in sorted(FIG_LEGENDS, key=lambda x: int(x[0].split('S')[1])):
        p = doc.add_paragraph()
        p.add_run(f'{name} - ').bold = True
        p.add_run(text)
    finalize(doc)
    doc.save(OUT)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
