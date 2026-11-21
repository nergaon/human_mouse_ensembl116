#!/usr/bin/env python3
"""Answers to Tal Shay's comments in Human_Mouse_Splicing_revised_tracked_TS2.docx (six-dataset version: three human and three mouse datasets per cell type).
Output: ensembl116/revision_6ds/manuscript/TS2_comment_responses.docx"""
from docx import Document
from docx_style_116 import apply_style, finalize

R = [
 ('General: design of this version', 'Each cell type is now compared between three human (H1-H3) and three mouse datasets (M1-M3), each with at least three replicates of the cell type (3-6 in mouse). The mouse datasets were assembled per cell type from 13 GEO series and one BioProject (GEO query in Methods 2.2, series published until 19 October 2026); GSE116177 (two replicates per cell type) is used only as a tissue control. Differentially spliced: LeafCutter ANOVA over the six datasets (adjusted P ≤ 0.05), then in at least one junction a human-mouse gap ≥ 0.1 and every human-mouse dataset difference > 1.5 × every within-species difference. Splicing conserved: ANOVA P > 0.1 and dataset PSIs within 0.05. The word "equivalent" was replaced by "splicing conserved" everywhere.'),
 ('Abstract: "How many of those that are testable in fibroblasts?"', 'Now stated: of the 470 events differentially spliced in at least one immune cell type and testable in fibroblasts, 45% are differentially spliced in fibroblasts.'),
 ('Introduction: "Did we lose a ref here?"', 'No; the submitted version had no reference on that sentence (refs 11-13 were on the next one). Added (11-13) to it.'),
 ('2.1: reasons for not using ExOrthist', 'Added: ExOrthist is restricted to protein-coding exons, whereas many alternative splicing events involve untranslated exons (e.g., alternative first/last exons).'),
 ('2.1: "Specify which" (arbitrary choice)', 'Specified the full tie-break order (canonical transcripts, then higher average identity, then coordinate order); coordinate order decided 4,794 of the 390,686 human boundaries (1.2%).'),
 ('2.1: "Are those results reported later?" (70/80%/RBH)', 'Yes: Results 3.1 (boundaries) and 3.3 (classification), Tables S1 and S6; pointer added.'),
 ('2.2: "All of these?" (adapters)', 'Yes, one cutadapt command removes all three (Nextera, TruSeq, Smart-seq2 ISPCR) from every sample of H3 and of the mouse datasets; rephrased.'),
 ('2.2: unrelated sentence after the moved controls paragraph', 'Removed; the controls are in Supplementary Note 1 (your text, with PRJNA797236 added).'),
 ('2.3: event types – supplementary figure? Do we use this classification?', 'The types are used for the counts in Results 3.1 and for two sensitivity analyses (splice-site correspondence; splice sites up to 6 nt apart); the text now says so. A schematic supplementary figure of the event types was not made yet.'),
 ('2.4: "label the clusters – with what?"', 'With gene names; added.'),
 ('2.4: Tables S4/S5 "Move to Excel"', 'Done: Table_S4.xlsx and Table_S5.xlsx; all supplementary tables are now Excel files.'),
 ('2.4: "±δ – what is that?" / "equivalent to anything that would not pass t-test?"', 'The equivalence test (90% CI within ±δ) was removed. In single comparisons (dataset pairs, fibroblasts, controls) an event is now splicing conserved if the LeafCutter P value is > 0.1 and max|ΔPSI| ≤ 0.05, i.e. not significant at a lenient threshold and a small difference; a non-significant event with a larger difference is inconclusive. Methods 2.4-2.5 rewritten.'),
 ('2.5: "unclear" (equivalent definition)', 'Replaced: "equivalent" is now "splicing conserved", defined as above for single comparisons, and for the human-mouse classification as ANOVA P > 0.1 and dataset PSIs within a range of 0.05 in every junction.'),
 ('2.5: Table S6 "not clear, let\'s discuss"', 'Kept; explanation of each sheet sent separately, to discuss.'),
 ('2.5: "How many events this removes? ... at least two junctions"', 'The PSI is computed from all reads of the event, so a sample has a PSI for all junctions or none; there is no case of one junction missing. Testability is now per event: a PSI in at least two samples of each of the six datasets (2,480-2,955 testable events per cell type).'),
 ('2.5: "Reconsider after LeafCutter all and Wilcoxon"', 'Done. Pooled LeafCutter and Wilcoxon alone also call many events in splits that mix the species (dataset effect). The classification is now: LeafCutter ANOVA over the six datasets (Dirichlet-multinomial model, all junctions, the six datasets as six groups, as in our 2023 paper), without the species labels. Splicing conserved: ANOVA P > 0.1 and range ≤ 0.05 in every junction. Differentially spliced: ANOVA adjusted P ≤ 0.05 and, in at least one junction, |mean human - mean mouse PSI| ≥ 0.1 and every human-mouse dataset |ΔPSI| > 1.5 × every within-species |ΔPSI|. The same rule applied to the nine splits that mix the species gives 1.8-4.1% as many events as the species split (Table S6; DSA tables). No confidence levels.'),
 ('2.6: "What? Why did you add this?" (coverage regression)', 'Removed: the coverage regressions (testability, classification and the coverage-adjusted enrichment) and their table were removed from the manuscript; Tables S8-S13 are now S7-S12.'),
 ('3.1: "unique – correct?"', 'Yes (exons are unique by coordinates across transcripts).'),
 ('3.1: reciprocal best gives more boundaries – strange', 'Explained in the text: removing non-reciprocal exon pairs removes competitors in the mutual-best pairing, so 1,478 boundaries that lost to a non-reciprocal rival become paired (390,630 shared).'),
 ('3.1: "Not shown" (97% constant offset)', 'It is in Table S1, sheet rank_consistency; pointer added.'),
 ('3.1: "This does not belong in this paragraph"', 'The testable-events sentence moved to 3.2.'),
 ('3.2: Figure 3A "not showing that. Maybe replace by a scatter of max ΔPSI?"', 'Not changed yet; Figure 3A still shows the cumulative distributions per pair. A scatter can be added if you prefer.'),
 ('3.2: "mean – correct?" / "replace the word confidence" / "max-min < 0.1?" / "By definition no?"', 'Method rewritten (ANOVA and the split rule), no confidence levels; range = max - min ≤ 0.1 stated; "as required by their definition" added.'),
 ('3.2/3.3: robustness "In what cell type? Robust in what way?"', 'Numbers per cell type and of unique events for 70%, 80% and reciprocal-best maps given.'),
 ('3.3: coverage effect "explain to me"', 'The coverage analysis was removed from the manuscript.'),
 ('3.3: SRSF7 "remove if no figure; specify in how many cell types"', 'Number of cell types given for each example (RPL22L1 3, PUF60 3, USP16 1); the 10 events DS in all five cell types are listed. SRSF7 is mentioned only in the Discussion (its A5SS is DS in four cell types).'),
 ('3.3: "Correct?" (single cell type) / "Both high and low?" (enrichment)', 'Rewritten with the new numbers (46% not testable, 54% inconclusive, no event splicing conserved). Only the pooled DS-only gene list was tested for enrichment (no confidence levels now); the added sentence about shared/each cell type lists was removed because those were not tested.'),
 ('3.4: gene examples "Which are what?" / different number / complicated / fraction of all testable / "Where was this shown?" / "Correct?"', 'Functional groups given; one consistent count (477 splicing conserved classifications, 4% of 13,425 testable); text simplified; Figure 3C/Table S6 cited for the minor-isoform result; "JSRs of the event" confirmed. A better visualization than Figure 3C is still open.'),
 ('3.4: "How was the adjustment done?" / "Pooling is per type?"', 'The coverage adjustment was removed. Pooling is per gene over all events and cell types; 41 genes with both classes are in neither list.'),
 ('3.5: move dataset details to methods; "Correct?"', 'Details are in Methods 2.2; removed from Results. Yes, in immune cells.'),
 ('Discussion: "How?" (orthology does not guarantee...)', 'Sentence rewritten to explain.'),
 ('Discussion: NAGNAG paragraph "weird now; characterize the list"', 'The 10 events DS in all five cell types are characterized by type and gene; NAGNAG kept for USP16/RPL22L1 with the test on all A3SS/A5SS.'),
 ('Discussion: "Requires ref"', 'Added Lek et al. 2016 (sequence constraint) and Brawand et al. 2011 (expression evolution).'),
 ('Discussion: "Of what?" / "Where is this shown?"', '"An example of an event that differs between species in immune cells but not in fibroblasts"; numbers from Table S11.'),
 ('Table 1: header one sentence, footnotes; adjusted p; n and % confusing', 'Caption is one sentence; definitions are numbered footnotes in the table; all classes are % with testable n rows in every section; "adjusted P".'),
 ('Figure 1: move splits to Note 2; may not pass reviewers; 5-junction event; bottom illustration', 'Splits removed from Figure 1; Figure 1B redrawn for the six datasets (mouse accessions per cell type) as the flow testable → splicing conserved (ANOVA P > 0.1, range ≤ 0.05) → differentially spliced (ANOVA, gap ≥ 0.1, 1.5× rule) → inconclusive; bottom illustration removed. Panel A still labels exons a/b (text now uses A/B).'),
 ('Figure 2: why not 100%, which pairs; unclear; "?"; remove gene symbols', 'The bars are now the mean over the dataset pairs, so each sums to 100% (per-pair values in Table S9); unclear sentence removed; "(Methods 2.5, Figure 1B)"; gene symbols removed from 2G.'),
 ('Figure 3: x-axis every 5-10%; colorbar; supp figure of all comparisons', 'Not done yet.'),
 ('Figure 6: numbers very different; testable %?', 'Legend now explains that the number above each bar is the denominator (events testable in fibroblasts and the cell type).'),
 ('"Supplementary figure legends???"', 'They are in Supplementary_Material.docx.'),
 ('Figure 5 (SRSF7 exon 7)', 'SRSF7 exon 7 skipping is splicing conserved only in NK cells, so Figure 5 is now the alternative 3\' splice site of EIF4B exon 7 (minor isoform <1% in both species; splicing conserved in four immune cell types, inconclusive in B cells, and in fibroblasts), a typical splicing conserved event.'),
]


def main():
    doc = Document()
    apply_style(doc)
    doc.add_heading('Responses to the comments in Human_Mouse_Splicing_revised_tracked_TS2', level=1)
    t = doc.add_table(rows=1, cols=2)
    t.style = 'Table Grid'
    t.rows[0].cells[0].text, t.rows[0].cells[1].text = 'Comment', 'Response'
    import table_numbers_116 as tn
    for c, r in R:
        cells = t.add_row().cells
        cells[0].text, cells[1].text = c, (r if 'Tables S8-S13 are now' in r else tn.renum(r))
    finalize(doc)
    out = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_6ds/manuscript/TS2_comment_responses.docx'
    doc.save(out)
    print('wrote', out)


if __name__ == '__main__':
    main()
