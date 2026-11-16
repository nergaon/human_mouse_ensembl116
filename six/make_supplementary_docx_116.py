#!/usr/bin/env python3
"""
Supplementary Material document of the revision (replaces the submitted Table_1.DOCX):
Tables S4, S5 and S10 as tables, the list of supplementary data sheets (Excel), and the legends of
Figures S1-S10.  Output: ensembl116/revision_6ds/supplementary/Supplementary_Material.docx
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
ED = B + 'ensembl116/events/id50_6ds_anova/'
OUT = B + 'ensembl116/revision_6ds/supplementary/Supplementary_Material.docx'
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
CELL_NAME = {'CD4T': 'CD4T', 'CD8T': 'CD8T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mono', 'Fibro': 'Fibroblasts'}
DATASETS = [('H1', 'Human', 'GSE60424'), ('H2', 'Human', 'GSE115736'), ('H3', 'Human', 'GSE118165'),
            ('M1', 'Mouse', ''), ('M2', 'Mouse', ''), ('M3', 'Mouse', ''),
            ('FH', 'Human', 'E-MTAB-5919'), ('FM', 'Mouse', 'E-MTAB-5919')]

DATA_SHEETS = [
    ('S1', 'Human-mouse exon orthology (Ensembl 116): orthologous exons and unique orthologous exon boundaries for '
           'orthologous exons defined by at least 50% (main), 70% or 80% identity, or by reciprocal best matches, and '
           'the exon-rank consistency of the boundaries used by the alternative splicing events.'),
    ('S2', 'Samples used: dataset, accession, cell type group, tissue, sorting markers, cell state, treatment, donor '
           'or mouse identity, sex, age, strain, read length and mapped reads.'),
    ('S3', 'Datasets of this study and candidate datasets from the dataset searches (Methods 2.2, published until 19 '
           'October 2026; control datasets, 19 October 2026), with a database link, whether each was '
           'used in the main and in the control analysis, and why not.'),
    ('S4', 'Datasets summary: dataset, accession, cell type, number of replicates and average mapped reads, and the '
           'LeafCutter -i and -g parameters.'),
    ('S5', 'Compared samples: one row per comparison of two groups of samples, with the reason for the comparison.'),
    ('S6', 'Sensitivity analyses: the human-mouse classification with each orthology definition; effect sizes on the '
           'events testable in all 15 dataset pairs; events with well-supported splice-site correspondence; '
           'classification by minor-isoform usage; alternative splice sites up to six nucleotides apart; and the '
           'status in the other cell types of events differentially spliced in one cell type.'),
    ('S7', 'Effect of coverage, number of junctions and number of events per gene on testability and classification '
           '(logistic regression; odds ratios and P values).'),
    ('S8', 'Orthologous junctions: each junction of the alternative splicing events with its orthologous mouse '
           'junction, the exon ranks at the junction ends, the event, the gene and the event type.'),
    ('S9', 'Single dataset pairs: for each pair of datasets and cell type, the testable events and the percentages of '
           'differentially spliced, splicing conserved and inconclusive events (Figure 2A-C).'),
    ('S10', 'Human-mouse classification of the alternative splicing events in each immune cell type: testable, '
            'differentially spliced, splicing conserved and inconclusive events.'),
    ('S11', 'Classification of each orthologous alternative splicing event in each immune cell type and in '
            'fibroblasts.'),
    ('S12', 'Functional enrichment: the significant terms of the genes with only splicing conserved events, and the terms that '
            'remain significant after adjusting for coverage.'),
    ('S13', 'Sex, age and tissue control comparisons (Supplementary Note 1): events, % differentially spliced, median '
            'max(|ΔPSI|) and PSI correlation in the control and in the human-mouse pairs, and the paired Wilcoxon test.'),
]

FIG_LEGENDS = [
    ('Figure S1', 'The status of an orthologous alternative splicing event in a single comparison of two datasets, '
                   'or of two sample groups, in one cell type. A. Differential splicing was analyzed in each immune cell '
                   'type between (i) each of the three human and each of the three mouse datasets (used only for effect sizes, '
                   'Figure 3), (ii) each pair of human datasets and (iii) each pair of mouse datasets, and (iv) between human '
                   'and mouse fibroblasts. B. The comparisons classified in this way (human-mouse events in immune cells '
                   'are classified only from all six datasets, Figure 1B), including the sex, age and tissue control comparisons, '
                   'and where their results are shown. C. The decision process: an event is not testable if fewer than the '
                   'minimal number of samples (g) in a group have at least 20 JSRs in the event; differentially '
                   'spliced if the LeafCutter Benjamini-Hochberg adjusted P value (p.adjust) is at most 0.05 and '
                   'max(|ΔPSI|) is at least 0.1; splicing conserved if the LeafCutter P value is above 0.1 and max(|ΔPSI|) is '
                   'at most 0.05; and inconclusive otherwise. The human-mouse classification of each event across all six '
                   'immune cell datasets is shown in Figure 1B.'),
    ('Figure S2', 'Classification of human-mouse orthologous alternative splicing events, all events. A-C. As in '
                  'Figure 2, but as percentages of all 6,784 orthologous alternative splicing events, including those that '
                  'are not testable in the comparison (fewer than the minimal number of samples with at least 20 JSRs in a '
                  'group; Figure S1C). D-F. As in Figure 2. G. Heatmap of the events that are differentially spliced (green) '
                  'or splicing conserved (blue) in exactly two immune cell types (events in three or more are shown in Figure 2G), '
                  'and their status in fibroblasts (Fb).'),
    ('Figure S3', 'Splicing differences of the classified events in every comparison. For each immune cell type, '
                  'the max(|ΔPSI|) of (A) the differentially spliced and (B) the splicing conserved events ('
                  'n, number of events) in each of the three human-human (blue), three mouse-mouse (green) and nine '
                  'human-mouse (orange) comparisons in which the event is testable. Dotted line, ΔPSI = 0.1; numbers '
                  'above the boxes, median max(|ΔPSI|). Boxes, interquartile range; line, median; whiskers, 5th-95th '
                  'percentiles.'),
    ('Figure S4', 'The orthologous skipped exon 5 event in PUF60 gene displays differential splicing between human '
                  'and mouse. A. The schematics of the orthologous skipped exon event. Exons are represented by '
                  'rectangles, splicing junctions by arcs. Exon numbers are those of the human Ensembl canonical transcript '
                  'and of the mouse transcript most similar to it that contains the orthologous junctions (Puf60-202, '
                  'ENSMUST00000100527; the mouse Ensembl canonical transcript lacks exons 2-5). B. Stacked bars display the number of junction spanning '
                  'reads that cover the alternative splicing junctions in each sample (top), and the same numbers '
                  'normalized to the total to display PSI values (bottom). Colors are as in A. M3, M2, M1 – mouse '
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
                  'mouse immune cells (CD8 T cells and monocytes), and is not differentially spliced between human and mouse fibroblasts. A. The schematics of '
                  'the orthologous skipped exon event. Exons are represented by rectangles, splicing junctions by '
                  'arcs. B. As in Figure S4B. Colors are as in A.'),
    ('Figure S7', 'Sex, age and tissue of origin cause much smaller splicing differences than species. Each '
                  'control comparison is between two groups of samples of the same dataset (same cell type, '
                  'laboratory and protocol), within one species (H, human; M, mouse), that differ in one factor: sex (H: GSE107011, male vs female donors; '
                  'M: ENCODE spleen B cells, male vs female mice), age (H: GSE143213, young vs old donors of naive CD4 T '
                  'cells; M: GSE184496, young vs old mice, naive CD4 T cells) and tissue (H: GSE133383, NK cells from '
                  'blood vs spleen of the same donors; M: GSE233347, NK cells from spleen vs liver; GSE116177, '
                  'monocytes from bone marrow vs blood; PRJNA797236, Ly6Chi monocytes from blood vs spleen), compared with the human-human, mouse-mouse and human-mouse '
                  'dataset pairs of the same cell type (mean over pairs), all tested as single comparisons (Figure S1C). '
                  'A. Percentage of differentially spliced events (the denominator includes the inconclusive and '
                  'splicing conserved events): for the control comparisons, of the events testable in the control comparison '
                  'and in all nine human-mouse pairs of the cell type; for the human-human, mouse-mouse and human-mouse '
                  'dataset pairs, of all events testable in each pair, as in Figure 2C (the human-mouse pairs are '
                  'shown for comparison; human-mouse events are classified from all six datasets, Table 1). B, C. On the events testable in the control comparison and in all nine '
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
    meta = pd.read_csv(B + 'ensembl116/revision_6ds/sample_metadata_separation.tsv', sep='\t')
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
    """Table S5: one row per comparison (group 1 vs group 2), taken from the LeafCutter comparisons that were run.
    Group names as in the group_name column of Table S2, which lists the samples of each group."""
    order_lab = ['H1', 'H2', 'H3', 'M1', 'M2', 'M3', 'FH', 'FM']
    order_cell = CELLS + ['Fibro']
    pairs = set()
    for line in open(ED + 'leafcutter/jobs.txt'):
        sample_set, pair, cell = line.strip().rstrip('/').split('/')[-3:]
        if sample_set == 'matched':
            continue
        pairs.add((cell, *sorted(pair.split('_'), key=order_lab.index)))

    def reason(cell, a, b):
        if cell == 'Fibro':
            return 'Human vs mouse in a non-immune cell type, both species from the same study'
        if a[0] == b[0]:
            sp = 'human' if a[0] == 'H' else 'mouse'
            who = 'donors' if sp == 'human' else 'mice'
            return f'Within-species baseline: differences between {sp} datasets (laboratory, protocol, {who})'
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
    # Supplementary Note 1 (Tal Shay, Supplementary_note_1.docx, accepted): sex, age and tissue controls
    from text_6ds_updates import SUPP_NOTE_1
    from revised_text_116 import NEW_REFERENCES
    doc.add_heading('Supplementary Note 1 – Sex, age and tissue effect', level=2)
    import table_numbers_116 as tn
    for para in SUPP_NOTE_1:
        doc.add_paragraph(tn.renum(para))
    doc.add_heading('Supplementary References', level=3)
    for k in ('dogra', 'encode', 'monaco'):
        doc.add_paragraph(NEW_REFERENCES[k])

    doc.add_heading('Supplementary Tables', level=2)
    doc.add_paragraph('All tables are provided as separate data sheets '
                      '(Excel files named Table_S<number>.xlsx).')
    import table_numbers_116 as tn   # published table numbers (the coverage table was removed)
    sheets = {tn.name(k): tn.renum(v) for k, v in DATA_SHEETS if tn.name(k)}

    def text_s3():
        add_table(doc, ['Dataset', 'Accession', 'Cell type', '#Replicates', 'Average mapped reads'], table_s1())
        doc.add_paragraph('-i (min_samples_per_intron) and -g (min_samples_per_group) of LeafCutter were set to 60% '
                          'of the smallest group, and at least 2: 3 for H2 vs M2 in CD4T, CD8T and NK, and 2 in all '
                          'other comparisons. Mapped reads: SAMtools flagstat mapped reads (alignments of paired reads '
                          'are counted per mate).')

    def text_s4():
        add_table(doc, ['Group 1', 'Group 2', 'Reason for the comparison'], compared_samples(),
                  widths=[3.5, 3.5, 9])
        doc.add_paragraph('Each row is one comparison of group 1 vs group 2. Group = dataset:cell type (datasets as in '
                          'Table S4); the samples of each group are listed in Table S2 (column group_name). '
                          'FH vs FM (fibroblasts) were compared '
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
                          'High + low confidence includes the high-confidence events. Also provided as Table_S10.xlsx.')

    # Tables S4 and S5 as Excel sheets (Tal Shay: move to Excel); S10 is written by build_supp_tables_116.py
    import build_supp_tables_116 as bst
    bst.write('Table_S4', 'Supplementary Table S4: Datasets summary',
              '-i (min_samples_per_intron) and -g (min_samples_per_group) of LeafCutter were set to 60% of the smallest '
              'group, and at least 2: 3 for H2 vs M2 in CD4T, CD8T and NK, and 2 in all other comparisons. Mapped '
              'reads: SAMtools flagstat mapped reads (alignments of paired reads are counted per mate).',
              {'datasets': pd.DataFrame(table_s1(), columns=['Dataset', 'Accession', 'Cell type', 'Replicates',
                                                             'Average mapped reads'])})
    bst.write('Table_S5', 'Supplementary Table S5: Compared samples',
              'Each row is one comparison of group 1 vs group 2. Group = dataset:cell type (datasets as in Table S4); '
              'the samples of each group are listed in Table S2 (column group_name). FH vs FM (fibroblasts) were '
              'compared as a single pair.',
              {'comparisons': pd.DataFrame(compared_samples(), columns=['Group 1', 'Group 2', 'Reason'])})
    text_tables = {}
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
        p.add_run(tn.renum(text))
    finalize(doc)
    doc.save(OUT)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
