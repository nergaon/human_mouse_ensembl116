#!/usr/bin/env python3
"""
Methods and legends for the candidate sample-correlation figures (figure_sample_correlation_116.py).
Also the events separated by each 3-vs-2 split of the datasets (figure_split_counts_116.py) and the status of the
events in each dataset pair (figure_pair_status_116.py).
Output: ensembl116/revision_separation/candidate_figures/Candidate_figures_methods.docx
"""
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt

from docx_style_116 import apply_style, finalize

D = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_separation/candidate_figures/'
OUTNAME = 'Candidate_figures_methods.docx'
SETS = [('all', 'All junctions'), ('one', 'One junction per event'), ('common', 'Junctions in all 104 samples')]


def n_info(mode):
    n = pd.read_csv(D + f'sample_psi_correlation_{mode}_n.tsv', sep='\t', index_col=0).to_numpy()
    off = n[~np.eye(len(n), dtype=bool)]
    return int(n.diagonal().max()), int(pd.Series(off).median()), int(off.min()), int(off.max())


def split_section(doc):
    df = pd.read_csv(D + 'split_counts.tsv', sep='\t')
    cells = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
    lab = {'CD4T': 'CD4 T', 'CD8T': 'CD8 T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo'}
    hm = df[df['human_mouse']].set_index('cell')
    null = df[~df['human_mouse']].assign(two=lambda x: x['split'].str.split(' \\| ').str[1])
    rng = lambda x: f'{int(x.min()):,}-{int(x.max()):,}'
    e = null.groupby('cell')[['events_gap_ge_0.1', 'events_gap_ge_0.2']].mean()
    fdr_l = e['events_gap_ge_0.1'] / hm['events_gap_ge_0.1']
    fdr_h = e['events_gap_ge_0.2'] / hm['events_gap_ge_0.2']
    big = {t: null[null['two'] == t] for t in ('H2 H3', 'H1 H2')}
    rest = null[~null['two'].isin(big)]
    ratio = hm['events_gap_ge_0.1'] / null.groupby('cell')['events_gap_ge_0.1'].max()

    doc.add_heading('2. Events separated by each split of the five datasets', level=2)
    doc.add_heading('Methods', level=3)
    for p in [
        'Splits. The five datasets of each cell type (H1, H2, H3, M1, M2) can be divided into a group of three and a '
        'group of two in ten ways (Figure 1B of the manuscript). One of them is the human-mouse split (H1 H2 H3 vs '
        'M1 M2); each of the other nine puts human and mouse datasets in the same group, so that what it separates '
        'reflects differences between datasets rather than between species.',
        'Separated events. The testable events and the dataset PSI (mean of the per-sample PSI over the samples '
        'with at least 20 junction-spanning reads in the event) were exactly those of the classification. For each '
        'split, an event was separated when, for at least one of its junctions, the PSI of every dataset of one '
        'group differed from the PSI of every dataset of the other group by at least 0.1 (gap ≥ 0.1; high and low '
        'confidence), or at least 0.2 (high confidence). For the human-mouse split, these are the differentially '
        'spliced events of Table 1.',
        'Estimated FDR. The expected number of events separated by dataset differences alone was the mean count '
        'over the nine mixed splits, and the FDR was this mean divided by the human-mouse count (as in Table 1).',
        'Software. Python 3.11 (pandas, matplotlib). Script: figure_split_counts_116.py, which uses the functions '
        'of classify_separation_116.py that produce Table 1.',
    ]:
        head, _, body = p.partition('. ')
        par = doc.add_paragraph()
        par.add_run(head + '. ').bold = True
        par.add_run(body)

    doc.add_heading('Results (brief)', level=3)
    doc.add_paragraph(
        f'Of the {rng(hm["testable_events"])} testable events in each cell type, the human-mouse split separated '
        f'{rng(hm["events_gap_ge_0.1"])} by a gap of at least 0.1 and {rng(hm["events_gap_ge_0.2"])} by at least '
        '0.2, more than any of the nine mixed splits in every cell type. The mixed splits separated a mean of '
        f'{e["events_gap_ge_0.1"].min():.0f}-{e["events_gap_ge_0.1"].max():.0f} events (gap ≥ 0.1; FDR '
        f'{fdr_l.min():.0%}-{fdr_l.max():.0%}) and {e["events_gap_ge_0.2"].min():.0f}-'
        f'{e["events_gap_ge_0.2"].max():.0f} events (gap ≥ 0.2; FDR {fdr_h.min():.0%}-{fdr_h.max():.0%}). The mixed '
        'splits were not equal: the two in which one human dataset is grouped with the two mouse datasets separated '
        f'many more events than the other seven: H1 with M1 and M2 (vs H2 H3), {rng(big["H2 H3"]["events_gap_ge_0.1"])} '
        f'events ({rng(big["H2 H3"]["events_gap_ge_0.2"])} with gap ≥ 0.2), and H3 with M1 and M2 (vs H1 H2), '
        f'{rng(big["H1 H2"]["events_gap_ge_0.1"])} ({rng(big["H1 H2"]["events_gap_ge_0.2"])}), compared with '
        f'{rng(rest["events_gap_ge_0.1"])} ({rng(rest["events_gap_ge_0.2"])}) for the other seven. Thus, H1 and H3 '
        'each differ from the other two human datasets in many events; these two splits contribute most of the '
        f'expected false positives, and the human-mouse split still separated {ratio.min():.1f}-{ratio.max():.1f} '
        'times as many events as the larger of them.')

    doc.add_paragraph().add_run('Table 2. Number of events separated by each split of the five datasets, gap ≥ 0.1 '
                                '(in brackets, gap ≥ 0.2).').bold = True
    splits = df['split'].drop_duplicates().tolist()
    t = doc.add_table(rows=1, cols=len(cells) + 1)
    t.style = 'Table Grid'
    for c, h in zip(t.rows[0].cells, ['Split (group of three vs group of two)'] + [lab[c] for c in cells]):
        c.text = h
        c.paragraphs[0].runs[0].bold = True
    idx = df.set_index(['split', 'cell'])
    for sp in splits:
        row = t.add_row().cells
        row[0].text = sp.replace(' | ', ' vs ') + (' (human vs mouse)' if sp.startswith('H1 H2 H3') else '')
        for k, c in enumerate(cells):
            v = idx.loc[(sp, c)]
            row[k + 1].text = f'{v["events_gap_ge_0.1"]} ({v["events_gap_ge_0.2"]})'
    extra = [('Mean of the 9 mixed splits',
              lambda c: f'{e.loc[c, "events_gap_ge_0.1"]:.1f} ({e.loc[c, "events_gap_ge_0.2"]:.1f})'),
             ('Estimated FDR', lambda c: f'{fdr_l[c]:.0%} ({fdr_h[c]:.0%})'),
             ('Testable events', lambda c: f'{hm.loc[c, "testable_events"]:,}')]
    for name, f in extra:
        row = t.add_row().cells
        row[0].text = name
        for k, c in enumerate(cells):
            row[k + 1].text = f(c)
    for r_ in t.rows:
        for c in r_.cells:
            for par in c.paragraphs:
                par.paragraph_format.space_after = Pt(0)
                par.paragraph_format.space_before = Pt(0)

    doc.add_heading('Figure legend', level=3)
    doc.add_picture(D + 'Figure_split_counts.png', width=Inches(6.5))
    par = doc.add_paragraph()
    par.add_run('Figure D. Number of events separated by each of the ten splits of the five datasets. ').bold = True
    par.add_run('Each group of bars on the x-axis is one of the ten ways to divide H1, H2, H3, M1 and M2 into a group '
                'of three and a group of two (Figure 1B); shaded, the human-mouse split (H1 H2 H3 vs M1 M2); the '
                'other nine mix the species. Bars, the number of events separated by the split in each cell type '
                '(colors; legend: number of testable events); an event is separated when, for at least one '
                'junction, the PSI of every dataset of one group differs from the PSI of every dataset of the other '
                'group by at least 0.1 (top) or 0.2 (bottom). Numbers above the human-mouse bars: events separated '
                '(the differentially spliced events of Table 1). Dashed line: mean over the nine mixed splits and '
                'the five cell types, the number expected from dataset differences alone.')


def pair_section(doc):
    d = pd.read_csv(D + 'pair_status.tsv', sep='\t')
    cells = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']
    lab = {'CD4T': 'CD4 T', 'CD8T': 'CD8 T', 'NveB': 'B', 'NK': 'NK', 'Mono': 'Mo'}
    pairs = ['H1_H2', 'H1_H3', 'H2_H3', 'M1_M2', 'H1_M1', 'H1_M2', 'H2_M1', 'H2_M2', 'H3_M1', 'H3_M2']
    k = d['comparison'].map(lambda c: 'X' if c[0] != c[3] else c[0])
    r = lambda col, kk: f'{d.loc[k == kk, col].min():.0f}-{d.loc[k == kk, col].max():.0f}%'
    rng = lambda sub, col: f'{sub[col].min():.0f}-{sub[col].max():.0f}%'
    h3 = d[d['comparison'].isin(['H3_M1', 'H3_M2'])]
    oth = d[(k == 'X') & ~d['comparison'].str.startswith('H3')]
    m12 = d[d['comparison'] == 'M1_M2']

    doc.add_heading('3. Status of the events in each dataset pair', level=2)
    doc.add_heading('Methods', level=3)
    for p in [
        'Dataset pairs. The five datasets of each cell type form ten pairs: three human-human (H1-H2, H1-H3, '
        'H2-H3), one mouse-mouse (M1-M2) and six human-mouse pairs. Figure 2A-C of the manuscript shows, for each '
        'kind of pair, the median over its pairs; here each pair is shown separately.',
        'Single-pair test. Each pair was compared by LeafCutter, as for Figure 2A-C (Methods 2.4-2.5, Figure S1C). '
        'An event was testable in a pair when enough samples of both datasets had at least 20 junction-spanning '
        'reads in the event (LeafCutter -g rule), so that the testable events differ between pairs. A testable event '
        'was differentially spliced when p.adjust ≤ 0.05 and max|ΔPSI| ≥ 0.1, equivalent when the 90% confidence '
        'interval of the ΔPSI of every junction was within ±0.1, and inconclusive otherwise. Percentages are of '
        'the testable events of the pair and cell type.',
        'Relation to the classification. The single-pair test is not used to classify events between species; '
        'the human-mouse pairs are a step before the classification from all five datasets (Table 1, Figure 2D), '
        'which requires the human-mouse difference to be present in every human-mouse pair and to exceed the '
        'differences between datasets of the same species.',
        'Software. Python 3.11 (pandas, matplotlib). Script: figure_pair_status_116.py.',
    ]:
        head, _, body = p.partition('. ')
        par = doc.add_paragraph()
        par.add_run(head + '. ').bold = True
        par.add_run(body)

    doc.add_heading('Results (brief)', level=3)
    doc.add_paragraph(
        f'Differentially spliced events were {r("pct_DS", "H")} of the testable events in the human-human pairs, '
        f'{r("pct_DS", "M")} in the mouse-mouse pair and {r("pct_DS", "X")} in the human-mouse pairs; equivalent '
        f'events were {r("pct_equivalent", "H")}, {r("pct_equivalent", "M")} and {r("pct_equivalent", "X")}. In '
        'every cell type, each of the six human-mouse pairs had a higher fraction of differentially spliced events '
        'than each of the human-human and mouse-mouse pairs. The pairs with H3 had the fewest differentially '
        f'spliced events among the human-mouse pairs (H3-M1 and H3-M2, {rng(h3, "pct_DS")}, compared with '
        f'{rng(oth, "pct_DS")} in the other four) and more inconclusive events ({rng(h3, "pct_inconclusive")} vs '
        f'{rng(oth, "pct_inconclusive")}), consistent with the lower coverage and higher variability of H3. The '
        f'mouse-mouse pair had {rng(m12, "pct_DS")} differentially spliced events, the lowest fraction in all '
        'cell types except monocytes, and, except in CD4 T cells, the largest number of testable events.')

    doc.add_paragraph().add_run('Table 3. Status of the testable events in each dataset pair: % differentially '
                                'spliced / % equivalent / % inconclusive (number of testable events).').bold = True
    t = doc.add_table(rows=1, cols=len(cells) + 1)
    t.style = 'Table Grid'
    for c, h in zip(t.rows[0].cells, ['Pair'] + [lab[c] for c in cells]):
        c.text = h
        c.paragraphs[0].runs[0].bold = True
    idx = d.set_index(['comparison', 'cell'])
    for pr in pairs:
        row = t.add_row().cells
        row[0].text = pr.replace('_', '-')
        for j, c in enumerate(cells):
            v = idx.loc[(pr, c)]
            row[j + 1].text = (f'{v["pct_DS"]:.0f} / {v["pct_equivalent"]:.0f} / {v["pct_inconclusive"]:.0f} '
                               f'({int(v["testable"]):,})')
    for r_ in t.rows:
        for c in r_.cells:
            for par in c.paragraphs:
                par.paragraph_format.space_after = Pt(0)
                par.paragraph_format.space_before = Pt(0)

    doc.add_heading('Figure legend', level=3)
    doc.add_picture(D + 'Figure_pair_status.png', width=Inches(6.5))
    par = doc.add_paragraph()
    par.add_run('Figure E. Status of the events in each of the ten dataset pairs. ').bold = True
    par.add_run('Each panel is one pair of datasets (title color: blue, human-human; green, mouse-mouse; orange, '
                'human-mouse), and each bar one cell type (T4, CD4 T cells; T8, CD8 T cells; B, B cells; NK, NK '
                'cells; Mo, monocytes). Bars show the percentage of the testable events of the pair that are '
                'differentially spliced (LeafCutter p.adjust ≤ 0.05 and max|ΔPSI| ≥ 0.1; green), equivalent (90% '
                'confidence interval of ΔPSI within ±0.1; blue) or inconclusive (gray); numbers in the bars are '
                'event counts (shown when the segment is at least 7%). Figure 2A-C of the manuscript shows the '
                'median of these percentages over the pairs of each kind.')


def main():
    s = pd.read_csv(D + 'sample_psi_correlation_summary.tsv', sep='\t', index_col=0)
    r = {m: pd.read_csv(D + f'sample_psi_correlation_{m}.tsv', sep='\t', index_col=0) for m, _ in SETS}
    total = {'all': 13049, 'one': 4569, 'common': 924}
    info = {m: n_info(m) for m, _ in SETS}
    lo = {m: r[m].to_numpy().min() for m, _ in SETS}

    doc = Document()
    apply_style(doc)
    doc.add_heading('Candidate figures', level=1)
    doc.add_heading('1. Correlation of junction PSI between samples', level=2)

    doc.add_heading('Methods', level=3)
    paras = [
        'Samples. All 104 immune-cell samples of the five datasets were used: H1 (GSE60424; 20 samples), H2 '
        '(GSE115736; 25), H3 (GSE118165; 20), M1 (24) and M2 (15), in five cell types (CD4 T, CD8 T, B and NK cells '
        'and monocytes). Fibroblasts were not included.',
        'Junction PSI per sample. The junction read counts of the orthologous alternative splicing events (the '
        'same junctions as in the main analysis, exon identity ≥50%) were used. In each sample, the PSI of a '
        'junction was its number of reads divided by the total number of reads of all junctions of its event '
        '(cluster), calculated only when the event had at least 20 junction-spanning reads in that sample, as in '
        'the classification of the events; otherwise the PSI was missing. Unlike the classification, PSI values '
        'were not averaged over the samples of a dataset: every sample is a row and a column of the matrix.',
        'Junction sets. The correlation was calculated with three sets of junctions. (1) All junctions: every '
        f'junction of every event ({total["all"]:,} junctions with a PSI in at least one sample). Each pair of '
        'samples was compared on the junctions that have a PSI in both samples (pairwise-complete), so that the '
        'junctions differ between pairs; in this set, the junctions of the same event are not independent (for an '
        'event with two junctions, the two PSI values sum to 1), and events with more junctions have more weight. '
        f'(2) One junction per event ({total["one"]:,} events): the junction with the lowest start coordinate '
        '(genomic, not transcript, order; if several junctions share that start, the first listed), again '
        'pairwise-complete. The junction was chosen by position only, independently of its PSI values or '
        'coverage, so that the choice does not favor high or low correlations; this removes the dependence between junctions of the same event and '
        'gives each event the same weight. (3) Junctions in all samples: the junctions with a PSI in all 104 samples '
        f'({total["common"]:,} junctions), so that all pairs are compared on the same junctions; these junctions are '
        'in highly expressed genes and include many events dominated by one isoform, with PSI values near 0 or 1, '
        'which increases the correlations.',
        'Correlation. For each pair of samples, the Pearson correlation coefficient (r) of the junction PSI values '
        'was calculated over the junctions of the set with a PSI in both samples (at least 100 junctions were '
        'required). The number of junctions per pair was median '
        f'{info["all"][1]:,} (range {info["all"][2]:,}-{info["all"][3]:,}) for all junctions, '
        f'{info["one"][1]:,} ({info["one"][2]:,}-{info["one"][3]:,}) for one junction per event, and '
        f'{total["common"]:,} for every pair with the junctions in all samples.',
        'Ordering of the samples. The samples were grouped by dataset and cell type (25 groups), and the replicates '
        'of each group were kept next to each other. The groups were ordered by average-linkage hierarchical '
        'clustering, with the distance between two groups defined as 1 minus the mean r over all pairs of their '
        'samples (dendrogram above the matrix). Within each group, the samples were ordered by average-linkage '
        'clustering of 1 - r of the group. White lines separate the groups. The bars above the matrix show the '
        'species, the dataset and the cell type of each sample. The color scale of each matrix starts at its lowest '
        f'r (0.54, 0.49 and 0.80 for the three sets), so that colors are not comparable between the three figures.',
        'Summary. For each junction set, the median r was calculated over the pairs of samples of each type: same '
        'dataset and same cell type (replicates), same dataset and different cell type, same species and different '
        'dataset, and human vs mouse; and for each pair of datasets (Table 1).',
        'Software. Python 3.11 with pandas (pairwise-complete Pearson correlation), SciPy (hierarchical clustering) '
        'and matplotlib. Script: figure_sample_correlation_116.py.',
    ]
    for p in paras:
        head, _, body = p.partition('. ')
        par = doc.add_paragraph()
        par.add_run(head + '. ').bold = True
        par.add_run(body)

    doc.add_heading('Results (brief)', level=3)
    doc.add_paragraph(
        'With all three junction sets, the samples were separated first by species and then by dataset; the cell '
        'types were ordered within each dataset. The median r between replicates of the same dataset and cell type '
        f'({s.loc["same dataset, same cell type", "one"]:.2f} with one junction per event) was only slightly higher '
        f'than between cell types of the same dataset ({s.loc["same dataset, other cell type", "one"]:.2f}), higher '
        f'than between datasets of the same species ({s.loc["same species, other dataset", "one"]:.2f}), which in '
        f'turn was higher than the human-mouse correlation ({s.loc["human vs mouse", "one"]:.2f}). The highest correlations were '
        f'within M1 ({s.loc["M1-M1", "one"]:.2f}), between M1 and M2 ({s.loc["M1-M2", "one"]:.2f}) and within H1 and '
        f'H2 ({s.loc["H1-H1", "one"]:.2f} and {s.loc["H2-H2", "one"]:.2f}; H1-H2 {s.loc["H1-H2", "one"]:.2f}). H3 '
        f'samples correlated poorly even with each other ({s.loc["H3-H3", "one"]:.2f}), consistent with the lower '
        'coverage of H3 (a median of 6,097 junctions with a PSI per sample, compared with 8,886-11,887 in the other '
        'datasets), but more with the other human datasets '
        f'({s.loc["H1-H3", "one"]:.2f} and {s.loc["H2-H3", "one"]:.2f}) than with the mouse datasets '
        f'({s.loc["H3-M2", "one"]:.2f}-{s.loc["H3-M1", "one"]:.2f}). The pattern was the same with the three junction sets '
        '(H1-M1 was as high as H1-H3 in all three); removing the redundant junctions lowered r by 0.01-0.06, and restricting to the '
        'junctions in all samples raised r in all pairs.')

    doc.add_paragraph().add_run('Table 1. Median Pearson r of junction PSI by type of sample pair and by pair of '
                                'datasets, for the three junction sets.').bold = True
    rows = [i for i in ['same dataset, same cell type', 'same dataset, other cell type', 'same species, other dataset',
                        'human vs mouse'] + [x for x in s.index if '-' in x]]
    t = doc.add_table(rows=1, cols=4)
    t.style = 'Table Grid'
    for c, h in zip(t.rows[0].cells, ['Sample pairs'] + [lab for _, lab in SETS]):
        c.text = h
        c.paragraphs[0].runs[0].bold = True
    for i in rows:
        cells = t.add_row().cells
        cells[0].text = i[0].upper() + i[1:] if not i[0].isupper() else i
        for k, (m, _) in enumerate(SETS):
            cells[k + 1].text = f'{s.loc[i, m]:.2f}'
    for row in t.rows:
        for c in row.cells:
            for par in c.paragraphs:
                par.paragraph_format.space_after = Pt(0)
                par.paragraph_format.space_before = Pt(0)

    doc.add_heading('Figure legends', level=3)
    common = ('Each row and column is one of the 104 immune-cell samples of the five datasets (H1, H2, H3, M1, M2); '
              'each cell is the Pearson correlation coefficient (r) between the junction PSI values of two samples '
              '(PSI calculated in each sample for events with at least 20 junction-spanning reads). Replicates of '
              'the same dataset and cell type are kept together and separated by white lines; the groups are ordered '
              'by hierarchical clustering of their mean correlation (dendrogram, average linkage, distance 1 - r), '
              'and the samples within a group by clustering of their own correlations. Bars above the matrix: '
              'species, dataset and cell type (T4, CD4 T cells; T8, CD8 T cells; B, B cells; NK, NK cells; Mo, '
              'monocytes). ')
    legends = {
        'all': f'All junctions of every event ({total["all"]:,} junctions); each pair of samples is compared on '
               f'the junctions with a PSI in both samples (median {info["all"][1]:,} junctions per pair).',
        'one': f'One junction per event (the first junction in genomic order; {total["one"]:,} events); each pair '
               f'of samples is compared on the events with a PSI in both samples (median {info["one"][1]:,} per pair).',
        'common': f'Only the {total["common"]:,} junctions with a PSI in all 104 samples, so that all pairs are '
                  'compared on the same junctions.',
    }
    for k, (m, lab) in enumerate(SETS):
        doc.add_picture(D + f'Figure_sample_correlation_{m}.png', width=Inches(6.0))
        par = doc.add_paragraph()
        par.add_run(f'Figure {"ABC"[k]}. Correlation of junction PSI between samples: {lab.lower()}. ').bold = True
        par.add_run(common + legends[m] + f' The color scale starts at the lowest r ({lo[m]:.2f}).')

    split_section(doc)
    pair_section(doc)
    finalize(doc)
    doc.save(D + OUTNAME)
    print('wrote', D + OUTNAME)


if __name__ == '__main__':
    main()
