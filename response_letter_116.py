#!/usr/bin/env python3
"""
Letter to the editor and point-by-point response to the reviewers (Frontiers in Immunology, 1982351).
Output: ensembl116/revision_files/Response_to_reviewers.docx and .md (the .md sections can be pasted into
the Frontiers review forum, one reviewer per thread).
"""
from docx import Document
from docx.shared import Pt, RGBColor

from docx_style_116 import apply_style, finalize

OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_files/Response_to_reviewers'

LETTER = [
    'Dear Editor,',
    'We thank you and the reviewers for the careful evaluation of our manuscript "Alternative splicing use ratio '
    'is not conserved between human and mouse immune cells" (manuscript ID 1982351). The reviewers raised two '
    'central concerns: that species differences are confounded with dataset, tissue and other biological factors, '
    'with only one dataset per species in the main comparison, and that "splicing unchanged" (no significant '
    'difference) was interpreted as conservation. We agree with both concerns, and we have substantially revised '
    'the study to address them. The main changes are:',
    '1. Additional independent datasets. The analysis now includes three human datasets (GSE60424, GSE115736 and '
    'GSE118165) and two mouse datasets (GSE116177, and a second mouse dataset assembled from four studies of '
    'sorted cells from unmanipulated or control mice, replacing the short-read GSE180020). Every human dataset is '
    'compared with every mouse dataset (six human-mouse pairs) and with each other, and the mouse datasets are '
    'compared with each other. The higher rate of human-mouse than within-species differential splicing is found '
    'in every human-mouse pair.',
    '2. Effect sizes independent of significance thresholds. We added distributions of |ΔPSI| and correlations of '
    'PSI for human-mouse, human-human and mouse-mouse comparisons on the same events (new Figure 3 and Figure S1), '
    'and analyses '
    'with alternative ΔPSI thresholds.',
    '3. Equivalence testing. "Splicing unchanged" was replaced by a three-way classification: differentially '
    'spliced, equivalent (90% confidence interval of ΔPSI within ±0.1, i.e. two one-sided tests) and '
    'inconclusive. Species divergence now requires a difference in two human-mouse dataset pairs that share no '
    'dataset and equivalence in no human-mouse pair; conservation requires demonstrated equivalence in two such '
    'pairs and no difference in any human-mouse pair. Events are assigned high confidence when, in addition, no '
    'within-species pair differs and an equivalent use ratio is demonstrated within each species. This '
    'revealed that an equivalent use ratio can be demonstrated mostly for events dominated by one isoform (the minor '
    'isoform accounts for <5% of the reads), which we now state explicitly.',
    '4. Robustness of the orthology map. The analysis was updated to Ensembl release 116, an error in the exon '
    'sequence extraction was corrected, and all analyses were repeated with 70% and 80% exon identity, with '
    'reciprocal best matches only, and for events with identical human-mouse splice-site geometry.',
    '5. Cell population matching. We provide a sample-level metadata table (Table S2) and repeated the analysis '
    'with the mouse populations restricted to naive T cells and blood monocytes.',
    '6. Sex, age and tissue. We added comparisons within each species that differ only in sex, age or tissue of '
    'origin (human and mouse). These caused far fewer and smaller splicing differences (0-6.6% of the events) '
    'than species (45-49% on the same events), and even smaller than the differences between datasets of the same '
    'species (new Results section 3.6, Figure S7, Table S13).',
    '7. Fibroblasts. We compared effect sizes on the same events and tested species x cell type interactions; the '
    'interpretation is now limited to supportive evidence, and we report that robust species differences are '
    'often shared with fibroblasts.',
    '8. Moderated conclusions. The title was changed to "Limited conservation of alternative splicing use ratio '
    'between human and mouse immune cells", the abstract was shortened and the denominators clarified, '
    'conclusions are restricted to sufficiently expressed and testable events, speculative evolutionary and '
    'functional statements were removed or marked as hypotheses, and the novelty of the approach is positioned '
    'relative to previous work.',
    'The central finding, that a large fraction of orthologous alternative splicing events differ in use ratio '
    'between human and mouse immune cells, more than between datasets of the same species, held in all of these '
    'analyses. A point-by-point response follows. Changes are marked in the tracked-changes version of the '
    'manuscript.',
    'Sincerely,',
    'Hadas Ner-Gaon, Naama Nudelman and Tal Shay',
]

# (reviewer, [(comment, response), ...])
RESPONSES = [
    ('Reviewer 1', [
        ('The major concerns are (1) potential confounding effects and (2) limited sample size. ... 1. Species '
         'differences are confounded with tissue source and other biological variables. The human-mouse '
         'comparisons do not isolate a species effect. ... human cells are mainly derived from peripheral blood, '
         'whereas mouse cells include spleen or other tissues. ... A tissue-source or age effect that is '
         'consistently aligned with species would remain classified as a high-confidence species difference.',
         'We thank the reviewer for this important point. We agree that a comparison of public human and mouse '
         'datasets cannot completely isolate a species effect, as tissue source, age and other factors are '
         'aligned with species in such data; this limitation is now stated explicitly in the Discussion. To '
         'reduce confounding as much as possible with public data, we made the following changes. (a) We added '
         'datasets and now compare three human and two mouse datasets that differ in laboratory, protocol, tissue '
         'source and mouse strain (Methods 2.2, Figure 1B, Table S2). The human-mouse rate of differential '
         'splicing is higher than the within-species rate in every one of the six human-mouse dataset pairs, in '
         'every cell type (Table 1, Table S8), so the difference is not attributable to a single pair of datasets. '
         '(b) We provide a sample-level metadata table (Table S2) with the tissue source, sorting markers, cell '
         'state, treatment, age, sex and strain of every sample. (c) We repeated the analysis with the mouse M1 '
         'samples restricted to naive T cells and blood monocytes, the populations most similar to the human '
         'samples; 39-49% (median) of the events remained differentially spliced between species, compared with '
         '48-52% in the full analysis (Table 1, Table S5). (d) The mouse datasets themselves span spleen, lymph '
         'nodes, bone marrow and peritoneal cavity, and the within-mouse comparison, which includes these tissue '
         'differences, shows a lower rate of differential splicing (14-27%) than any human-mouse comparison. We '
         'agree that factors that are aligned with species in all datasets, such as housing conditions or age, '
         'cannot be excluded, and we now discuss them as limitations (Discussion, limitations paragraph). (e) Most '
         'importantly, we directly tested whether sex, age or tissue of origin can cause splicing differences of '
         'the observed magnitude, using comparisons within each species that differ only in one of these factors, '
         'with the same cell type, laboratory and protocol: male vs female (human: five cell types, GSE107011; '
         'mouse: B cells, ENCODE), young vs old (human and mouse naive CD4 T cells, GSE143213 and GSE184496) and '
         'two tissues (human NK cells from blood and spleen of the same donors, GSE133383; mouse NK cells from '
         'spleen and liver, GSE233347; mouse monocytes from bone marrow and blood, GSE116177). On the same events, '
         'only 0-0.1%, 0-2.3% and 0-6.6% of the events were differentially spliced by sex, age and tissue, '
         'respectively, compared with 45-49% between human and mouse; the effect sizes were also much smaller '
         '(median max(|ΔPSI|) 0.018-0.044 vs 0.097-0.135) and the PSI correlations higher (0.98-0.995 vs '
         '0.79-0.87), and the human-mouse differences were larger than those of every control comparison (paired '
         'Wilcoxon test, P < 0.001; new Results section 3.6, Figure S7, Table S13). These factors therefore cannot '
         'account for the human-mouse differences.'),
        ('2. Insufficient dataset-level and biological replication. The central human-mouse comparison includes '
         'only one dataset per species ... several groups contain only two biological replicates ... the reported '
         '62-67% differential-splicing rate and the classification of nonsignificant events as "splicing '
         'unchanged" cannot be considered robust or generalized without validation using additional independent '
         'datasets and biological replicates.',
         'We agree. (a) The central comparison now includes three human and two mouse datasets (six human-mouse '
         'pairs), and the main results are reported as the median and range over these pairs (Table 1). The '
         'second mouse dataset of the original submission (GSE180020) was replaced, as its short reads (2x38 '
         'nucleotides) left many junctions uncovered; the new datasets have reads of 75-150 nucleotides and 2-4 '
         'replicates per cell type. (b) Species divergence now requires differential splicing in two human-mouse '
         'pairs that share no dataset (independent replication) and equivalence in no human-mouse pair, and high '
         'confidence also requires no differential splicing and an equivalent use ratio within each species. (c) '
         'We no longer '
         'classify nonsignificant events as "unchanged": we test for equivalence (90% confidence interval of '
         'ΔPSI within ±0.1, equivalent to two one-sided tests), so that groups with few replicates, which have '
         'wide confidence intervals, are classified as inconclusive rather than unchanged (Methods 2.4-2.5). '
         'With this definition, an equivalent use ratio could be shown for a median of only 11-18% of the events '
         'in the human-mouse comparisons, and the conclusions about conservation were revised accordingly (Results 3.4). '
         '(d) The human-mouse rate of differential splicing in the revised analysis is 48-52% (median over the six '
         'pairs) and is higher than the within-species rate (14-34%) in every pair.'),
        ('Statistical methods (No) and documentation for replication (No).',
         'We revised the statistical framework as described above (equivalence testing, replication across '
         'dataset pairs, effect-size analyses, coverage models) and expanded the Methods. All analysis scripts, '
         'including the exact versions of the annotation (Ensembl release 116, pinned to the release-116 files) '
         'and of the software, are available in the GitHub repository, and all intermediate results are provided '
         'in the Supplementary Tables (S4-S13).'),
    ]),
    ('Reviewer 2', [
        ('1. Species effects are strongly confounded with dataset and sample-source effects. ... The authors '
         'should more quantitatively separate the species effect from within-species dataset variability. For '
         'example, distributions of absolute ΔPSI values in human-mouse, human-human, and mouse-mouse comparisons '
         'should be directly compared, rather than relying mainly on the proportion of events passing '
         'significance thresholds.',
         'We thank the reviewer for this suggestion, which we implemented. In the new Figure 3, we compare the '
         'distributions of max(|ΔPSI|) between human-mouse, human-human and mouse-mouse comparisons on the same '
         'events (the 1,752-1,982 events per cell type that are comparable in all ten comparisons). In every cell '
         'type, human-mouse differences are larger than within-species differences (paired Wilcoxon test, '
         'P < 10^-94; Figure 3A). The distributions of max(|ΔPSI|) and of the 90% CI of ΔPSI in each of the ten '
         'comparisons are shown separately for each cell type in the new Figure S1. The human-mouse differences '
         'are larger on average, although comparisons that include the datasets with fewer replicates or '
         'Smart-seq2 libraries (H3 and M2) show larger differences and wider confidence intervals also within '
         'species. The rate of differential splicing is higher in every human-mouse dataset pair '
         'than in every within-species pair (Table S8). In addition, we added more datasets so that dataset variability is estimated from '
         'three human-human pairs and one mouse-mouse pair (Figure 1B), and we repeated the analysis with matched '
         'cell populations (Table S5), and we added control comparisons within each species that differ only in sex, '
         'age or tissue of origin, which showed far smaller differences than species (0-6.6% vs 45-49% of the events '
         'differentially spliced; Results 3.6, Figure S7, Table S13). We also discuss the remaining confounders (sequencing depth, library '
         'preparation, tissue source, age, sex, genetic background) in the Discussion.'),
        ('2. The headline conclusion based on 62-67% differential splicing is difficult to reconcile with the '
         'high-confidence analysis. ... The manuscript should therefore clarify which analysis constitutes the '
         'primary evidence supporting the biological conclusion. If low-confidence events are substantially '
         'affected by technical or dataset variability, the title, abstract, and overall interpretation should be '
         'moderated accordingly.',
         'We agree that this needed clarification. In the revised manuscript, the primary evidence for species '
         'divergence is the event classification, which requires differential splicing in two human-mouse dataset '
         'pairs that share no dataset (2,424 events, 1,164-1,389 per cell type), and, for high confidence, no '
         'difference and an equivalent use ratio within each species (120 events, 6-51 per cell type; Results 3.2, '
         'Table S9). For the differentially spliced events, the differences between datasets of the same species are '
         'smaller than the human-mouse differences (median max(|ΔPSI|) 0.06-0.16 vs 0.17-0.29; new Figure S3). The '
         'small number of high confidence events reflects the difficulty of demonstrating equivalence within species '
         'for events in which both isoforms are used and the frequent differential splicing between datasets of the '
         'same species, and we state this explicitly. The per-comparison rates (48-52% between species vs 14-34% within species) and the '
         'threshold-free effect sizes (Figure 3) are presented as supporting evidence that the species difference '
         'exceeds dataset variability on the whole set of testable events. We moderated the title ("Limited '
         'conservation of alternative splicing use ratio between human and mouse immune cells"), and the abstract '
         'now reports both the per-comparison rates and the number of differentially spliced events.'),
        ('3. "Splicing unchanged" should not automatically be interpreted as evidence of conservation. ... The '
         'authors should distinguish more clearly between "no evidence of differential splicing" and "evidence of '
         'conserved splicing." An equivalence-based framework, confidence intervals for ΔPSI, or at minimum a more '
         'cautious terminology would improve the statistical interpretation.',
         'We fully agree and adopted an equivalence-based framework. For each event and comparison we estimate '
         'the 90% confidence interval of the ΔPSI of every junction (Welch), and an event is equivalent only if all '
         'intervals lie within ±0.1 (equivalent to two one-sided tests at α = 0.05). Events that are neither '
         'differentially spliced nor equivalent are inconclusive (Methods 2.4-2.5, Figure 1C). Conservation '
         'requires equivalence in two human-mouse pairs that share no dataset and no differential splicing in any '
         'human-mouse pair, and high confidence also no difference and equivalence within each species. This had '
         'a substantial effect: equivalence could be demonstrated for a median of only 11-18% of the events '
         'between species, and 92% of the equivalent events are dominated by one isoform (the minor isoform accounts for <5% of the '
         'PSI; Figure 3D). We revised Results 3.4 and the Discussion accordingly, and we report results for '
         'margins of ±0.025 and ±0.05 as well (Table S8, S5).'),
        ('4. The definition of orthologous exons requires stronger validation and sensitivity analysis. ... The '
         'authors should provide sensitivity analyses using more stringent thresholds, for example 70% or 80% '
         'identity, and evaluate whether the principal conclusions remain stable. Incorporating additional '
         'information such as conserved genomic order, synteny, exon rank, or reciprocal matching would also '
         'strengthen the orthology definition.',
         'We repeated the entire analysis, from the orthology map to the classification, with orthologous exons '
         'defined by at least 70% or 80% identity, and by reciprocal best matches only (Table S1, S5). The '
         'principal conclusions are unchanged: the range of human-mouse differential splicing is 38-64% of the '
         'comparable events per pair at all stringencies, and 112, 88 and 122 events are differentially spliced with '
         'high confidence, compared with 120 at 50% identity. Synteny is incorporated through the use of high '
         'confidence one-to-one orthologous genes only (Ensembl orthology confidence, which includes gene order '
         'conservation and whole genome alignment), junctions are defined only between orthologous exon '
         'boundaries, and we added a test of splice-site correspondence: in 3,065 events the '
         'distance between the alternative splice sites, or the length of the skipped exon, is identical in human '
         'and mouse, and the results in this subset are the same (Table S5). Exon order is also supported by the '
         'data: a junction is analyzed only if both of its exons are orthologous and the junction connecting them is '
         'observed in both species, that is, the two exons are sequential in at least one expressed transcript of '
         'both human and mouse. Although exon rank is not used to define orthology, 85% of the orthologous exon '
         'boundaries used by the events have the same exon rank in human and mouse, and 97% have a constant rank '
         'offset within the gene (Table S1). While updating the pipeline to '
         'Ensembl 116 we also found and corrected an error in the extraction of exon sequences (line breaks of '
         'the sequence files were included in some aligned sequences, and exons with untranslated regions on both '
         'sides were trimmed on one side only).'),
        ('5. The fibroblast comparison does not fully establish immune-cell-specific divergence. ... I suggest '
         'interpreting the fibroblast results as supportive evidence rather than definitive evidence for '
         'immune-cell-specific divergence.',
         'We agree and revised this section (Results 3.5, Discussion). With the revised analysis, 3,801 of the '
         '4,570 events (83%) are comparable in fibroblasts. The fibroblasts are compared with the differentially spliced '
         'and equivalent events of each immune cell type, which are defined across all six human-mouse pairs: 32-38% '
         'of the differentially spliced events are differentially spliced also in fibroblasts, and 89-93% of the '
         'equivalent events are equivalent in fibroblasts (Figure 6A, B). On the same events, human-mouse differences are '
         'smaller in fibroblasts than in immune cells (median 0.03 vs 0.10-0.13), but the events that are '
         'differentially spliced with high confidence in immune cells are often differentially spliced also in '
         'fibroblasts (43%). We therefore no longer claim that most differences are immune-specific; instead we '
         'conclude that robust differences are often shared between cell types, while others may reflect '
         'immune-specific regulation, and we describe the fibroblast results as supportive rather than definitive '
         'evidence, noting the single dataset per species and the cultured state of the fibroblasts. The section '
         'title was changed to "Human-mouse splicing differences are larger in immune cells than in fibroblasts".'),
        ('6. Additional effect-size-based analyses would substantially strengthen the manuscript. ... the authors '
         'could report correlations between human and mouse PSI values across orthologous events and compare these '
         'with human-human and mouse-mouse correlations. Similarly, density plots, cumulative distributions, or '
         'violin plots of |ΔPSI| would allow readers to assess the magnitude of cross-species divergence '
         'independently of statistical significance.',
         'Thank you; these analyses are now in the new Figure 3: cumulative distributions of max(|ΔPSI|) for every '
         'dataset pair (Figure 3A), PSI correlations for every dataset pair and cell type (Figure 3B; human-mouse '
         'r = 0.76-0.89 vs 0.88-0.97 within species on the same events) and PSI scatter plots (Figure 3C). Box plots '
         'of max(|ΔPSI|) and of the 90% CI of ΔPSI for each comparison and cell type are in the new Figure S1. The '
         'values for all comparisons are in Tables S8 and S5.'),
        ('7. The scope of the conclusions should reflect the limited number of testable events. ... The '
         'conclusions should explicitly refer to sufficiently expressed and testable orthologous alternative '
         'splicing events in the analyzed datasets.',
         'We agree and revised the abstract and Discussion to refer explicitly to sufficiently expressed and '
         'testable orthologous alternative splicing events (4,570 events, 3,722 of them testable in two '
         'independent human-mouse pairs in at least one cell type). We also show that coverage affects testability and '
         'classification (Table S6).'),
        ('Minor 1. The abstract appears repetitive and should be shortened.',
         'The abstract was rewritten and shortened, and each number is reported once.'),
        ('Minor 2. The manuscript should consistently distinguish "alternative splicing event," "junction," '
         '"PSI," and "splicing ratio."',
         'We defined the terms in the Methods (2.3-2.4) and used them consistently in the revised sections: an '
         'alternative splicing event is a LeafCutter cluster of junctions; PSI is the fraction of the junction '
         'reads of an event that support a junction; and the use ratio of an event is described by the PSI of its '
         'junctions.'),
        ('Minor 3. The rationale for the selected ΔPSI thresholds should be explained and, ideally, supported by '
         'sensitivity analyses using alternative effect-size cutoffs.',
         'The rationale is now given in Methods 2.5, and we repeated the analysis with ΔPSI thresholds of 0.05 and '
         '0.2 (and equivalence margins of 0.025 and 0.05). The human-mouse rate is about twice the within-species '
         'rate at every threshold (0.05: 60-64% vs 18-43%; 0.2: 27-33% vs 5-18%; Results 3.2, Table S5).'),
        ('Minor 4. The functional enrichment analysis collapses event-level classifications to the gene level. ... '
         'the implications of this simplification should be discussed.',
         'Only 35 genes contain both a differentially spliced and an equivalent event with high confidence, so pooling '
         'events per gene has little effect on the gene lists; this is now stated in Results 3.4, and the '
         'event-level classification of every event is provided in Table S10.'),
        ('Minor 5. The enrichment of conserved events in ribosomal and translation-related genes may partially '
         'reflect expression and junction-coverage biases. The authors should consider whether gene expression '
         'level, number of testable junctions, or number of alternative splicing events per gene influences the '
         'probability of classification.',
         'We tested this directly. In a logistic regression, higher event coverage increases the probability that '
         'an event is testable and classified as equivalent, and decreases the probability that it is classified as '
         'differentially spliced; events with more junctions are less likely to be classified as equivalent (Table S6). '
         'Several RNA processing and respiration terms remained associated with equivalent events after '
         'adjusting for coverage (for example rRNA processing, adjusted odds ratio 2.85, P = 0.006; ribonucleoprotein complex biogenesis, '
         '1.91, P = 0.005; oxidative phosphorylation, 2.26, P = 0.007; Table S11). We also show that equivalent '
         'events are mostly '
         'events dominated by one isoform (the minor isoform accounts for <5% of the reads), and discuss that equivalence is easier to demonstrate for such events '
         '(Results 3.4, Discussion).'),
        ('Minor 6. Several grammatical and typographical issues remain throughout the manuscript ... formatting '
         'artifacts in the reported percentages and several spelling errors in the figure legends.',
         'We corrected the typographical errors and formatting artifacts that we identified (for example '
         '"PairwieseAligner", "than" for "then", "differntially" in the figure legends, and the split '
         'percentages in Table 1), and rewrote the figure legends.'),
    ]),
    ('Reviewer 4', [
        ('1. Definition of "splicing unchanged". ... The authors should distinguish evidence of differential '
         'usage, evidence of sufficiently similar usage, and inconclusive results, using an appropriately '
         'justified equivalence criterion or uncertainty estimates. The reported conservation rates and downstream '
         'analyses should then be reassessed using the revised classification.',
         'We agree and implemented exactly this three-way classification: differentially spliced, equivalent '
         '(90% confidence interval of ΔPSI within ±0.1 for every junction, equivalent to two one-sided tests; the '
         'margin equals the ΔPSI threshold for differential splicing, and results for ±0.025 and ±0.05 are '
         'reported) and inconclusive (Methods 2.4-2.5, Figure 1C). All downstream analyses were repeated with the '
         'revised classification: the differentially spliced and equivalent events and their confidence levels (Figure 2, Table S9), the cell type sharing '
         '(Figure 2E-G), the enrichment (Table S11) and the fibroblast comparison (Figure 6). Equivalence could be '
         'shown for a median of only 11-18% of the events between species, and mostly for events dominated by one '
         'isoform (the minor isoform accounts for <5% of the reads) (Figure 3D, Results 3.4).'),
        ('2. Comparability of cell populations across datasets. ... Please provide a sample-level metadata table '
         'specifying the actual samples included, their biological sources, sorting criteria, cell states, '
         'treatment conditions, and biological replicate identities. Analyses restricted to the most closely '
         'matched populations would help determine whether the principal findings persist after reducing these '
         'sources of confounding.',
         'Thank you. Table S2 now lists every sample with its accessions, biological source, sorting markers, '
         'cell state, treatment condition, donor or mouse identity, sex, age, strain, library type, read length '
         'and mapping statistics. As noted by the reviewer, not all populations are naive: the human CD4 and CD8 T '
         'cells of H1 and H2 are total populations, and the M1 groups include memory and effector T cells and '
         'inflammatory and bone marrow monocytes. We corrected the description in the Methods and repeated the '
         'analysis with M1 restricted to naive CD4 T cells, naive CD8 T cells and blood monocytes; the human-mouse '
         'rate of differential splicing remained 39-49% (median), compared with 48-52% (Table 1, Table S5). In '
         'addition, comparisons within each species that differ only in sex, age or tissue of origin (for example '
         'NK cells from blood and spleen of the same human donors) showed far fewer and smaller differences than '
         'species (Results 3.6, Figure S7, Table S13).'),
        ('3. Validation of exon and splice-junction orthology assignments. ... Please justify the 50% nucleotide '
         'identity threshold and clarify how it relates to the criteria used in the cited orthology methods. '
         'Benchmarking against independently supported mappings, reporting ambiguous assignments, and testing '
         'alternative thresholds would strengthen confidence in the analysis. The principal findings should also '
         'be evaluated in a subset of events with particularly well-supported splice-site correspondence.',
         'We clarified the definition (at least 50% of each of the two exons aligned to identical nucleotides) '
         'and its relation to ExOrthist: our nucleotide criterion is less strict than ExOrthist, which aligns exons '
         'at the protein level and requires conserved intron positions and neighboring exons; in our framework, the '
         'exon context is instead required through the RNA-seq data, as a junction is analyzed only if its two '
         'orthologous exons are sequential in at least one expressed transcript of both species (Methods 2.1). '
         'Consistent with this, 85% of the orthologous exon boundaries used by the events have the same exon rank '
         'in human and mouse, and 97% a constant rank offset within the gene (Table S1). In addition, a human and a mouse '
         'exon boundary are now paired only if each is the best-scoring candidate of the other (Methods 2.1); this '
         'prevents a mouse junction from being assigned to a non-orthologous human splice site when its best '
         'human partner is already paired with another mouse boundary. We also repeated all analyses with 70% and 80% identity and with reciprocal best matches only, which removes '
         'ambiguous many-to-one assignments (Table S1). The results are stable (112, 88 and 122 high confidence '
         'differentially spliced events vs 120; human-mouse differential splicing 38-64% per pair at all '
         'stringencies; Table S5). We also evaluated the principal findings in the 3,065 events with '
         'well-supported splice-site correspondence, in which the distance between the alternative splice sites '
         '(alternative 3\' or 5\' splice sites) or the skipped exon length is identical in human and mouse: 47-50% '
         'of these events are differentially spliced between species vs 24-28% within species, the same as in '
         'all events (Table S5). The orthologous mouse coordinates of every junction are provided in Tables S7 '
         'and S6.'),
        ('4. Evidence supporting immune specificity. ... statistical significance in one cell type but not another '
         'does not demonstrate a significant difference between their species effects. The authors should compare '
         'effect sizes for the same assessable events and, where the design permits, test species-by-cell-type '
         'interactions. Otherwise, the interpretation should be restricted to differences observed relative to the '
         'fibroblast samples examined.',
         'We performed both analyses on the same events (Results 3.5, Figure 6E-F, Table S12). Human-mouse '
         'differences are smaller in fibroblasts than in immune cells (median 0.03 vs 0.10-0.13, P < 10^-185), '
         'and the human-mouse ΔPSI in immune cells and fibroblasts is only moderately correlated (r = 0.46-0.52). '
         'A test of the species x cell type interaction (difference between the human-mouse ΔPSI in the immune '
         'cell type and in fibroblasts) is significant for 27-34% of the events and 43-55% of the events '
         'differentially spliced in the immune cell type. At the same time, events that are differentially spliced with '
         'high confidence are often differentially spliced also in fibroblasts. We therefore restricted the interpretation to '
         'differences relative to the fibroblast samples examined, and removed the claim that most differences '
         'are immune-specific.'),
        ('5. Biological and evolutionary interpretations require stronger support. ... The discussion of '
         'cancer-related evolutionary selection is particularly speculative and should be removed or clearly '
         'identified as an untested hypothesis. Independent validation of representative splicing events would '
         'strengthen the study, while functional or mechanistic claims should either receive appropriate '
         'experimental support or be revised to reflect their exploratory nature.',
         'We removed the discussion of cancer-related evolutionary selection, and the statements on selection for '
         'an isoform ratio and on effects on protein function are now presented as hypotheses that remain to be '
         'tested (Discussion). We also no longer state that cell type dependence excludes cis-regulatory '
         'differences; instead, we discuss that robust differences shared with fibroblasts may reflect '
         'cis-regulatory or splicing code changes. Regarding validation, the requirement that a species difference '
         'be reproduced in two human-mouse dataset pairs that share no dataset, generated by different '
         'laboratories and protocols, and be absent from every within-species comparison, provides independent '
         'replication of each differentially spliced event in the RNA-seq data; experimental validation of individual events is beyond the scope of this '
         'computational study, and we state this limitation.'),
        ('6. Positioning of the methodological novelty. ... The authors should specify which components of their '
         'mapping strategy represent methodological advances and which constitute the integration or application '
         'of existing approaches.',
         'We revised the last paragraph of the Introduction: exon-level orthology, annotation-free local splicing '
         'quantification and cross-species splicing comparisons have precedents, and our contribution is their '
         'integration into an exon-boundary orthology map that converts mouse junctions to human coordinates, so '
         'that the same event is tested directly in both species, together with its application to sorted immune '
         'cell populations profiled in several independent datasets per species, which we provide as a resource '
         '(Tables S7, S10).'),
        ('Minor 1. Clarify the denominators in the abstract and main results ... Please also reconcile the "67%" '
         'statement in the submission abstract with the "62-67%" range in the manuscript abstract and remove the '
         'duplicated text in the submission abstract.',
         'The abstract was rewritten with explicit denominators: 4,570 orthologous events in total, and per cell '
         'type the percentage of the events testable in each human-mouse comparison. The submission abstract was '
         'replaced by the revised manuscript abstract, without duplicated text.'),
        ('Minor 2. Correct the denominator accompanying Figure 5C. The statement that 141-163 events constitute '
         '81-85% of 269 events is numerically inconsistent.',
         'Thank you. The fibroblast results (now Figure 6) are reported with the correct conditional denominators: '
         'for each panel, the number of events comparable in both the immune cell type and fibroblasts is stated '
         'in the text and above each bar.'),
        ('Minor 3. Correct figure references, panel labels, and statistical notation ... Table 1 should clearly '
         'distinguish the human-human and mouse-mouse comparisons and specify whether the reported significance '
         'threshold refers to adjusted P-values.',
         'We corrected the figure references and panel labels; the gene figures now have a single panel B '
         '(counts and PSI), described accordingly in the legends. Table 1 now has separate blocks for the '
         'human-mouse, matched human-mouse, human-human and mouse-mouse comparisons, and the caption states that p '
         'is the Benjamini-Hochberg adjusted P value.'),
    ]),
]


def main():
    doc = Document()
    apply_style(doc)
    doc.add_heading('Response to the reviewers – manuscript 1982351', level=1)
    md = ['# Response to the reviewers – manuscript 1982351', '']
    doc.add_heading('Letter to the editor', level=2)
    md += ['## Letter to the editor', '']
    for p in LETTER:
        doc.add_paragraph(p)
        md += [p, '']
    for reviewer, items in RESPONSES:
        doc.add_heading(reviewer, level=2)
        md += [f'## {reviewer}', '']
        for comment, response in items:
            p = doc.add_paragraph()
            r = p.add_run('Comment: ' + comment)
            r.italic = True
            r.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
            p = doc.add_paragraph()
            p.add_run('Response: ').bold = True
            p.add_run(response)
            md += [f'*Comment: {comment}*', '', f'**Response:** {response}', '']
    finalize(doc)
    doc.save(OUT + '.docx')
    open(OUT + '.md', 'w').write('\n'.join(md))
    print('wrote', OUT + '.docx')


if __name__ == '__main__':
    main()
