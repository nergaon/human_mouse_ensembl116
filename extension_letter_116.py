#!/usr/bin/env python3
"""Letter to the editor requesting an extension of the revision deadline (manuscript 1982351).
Output: ensembl116/revision_files/Letter_to_editor_extension.docx and .md"""
from docx import Document

from docx_style_116 import apply_style, finalize
from docx.shared import Pt

OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_files/Letter_to_editor_extension'

PARAGRAPHS = [
    ('h', 'Request for an extension of the revision deadline – manuscript 1982351'),
    ('p', 'Dear Editor,'),
    ('p', 'Thank you for the opportunity to revise our manuscript "Alternative splicing use ratio is not conserved '
          'between human and mouse immune cells" (manuscript ID 1982351), and the reviewers for their thorough and '
          'constructive comments. We intend to address all of the points raised by the reviewers. Several of them '
          'require new data and new analyses, and we are writing to ask for an extension of the revision deadline.'),
    ('p', 'The reviewers\' main concerns were that the species difference may be confounded with dataset, tissue, '
          'age, sex and other biological factors, with only one dataset per species in the main comparison, and that '
          'the absence of a significant difference was interpreted as conservation. To address these concerns, we '
          'have already completed the following:'),
    ('b', 'Add independent datasets, so that the analysis includes three human (adding GSE118165) and two mouse '
          'RNA-seq datasets of the same five immune cell types. The mouse dataset GSE180020 will be replaced by a '
          'dataset assembled from GSE287527 (B cells), GSE261448 (NK cells and monocytes), ENCODE (ENCSR968ZOA and '
          'ENCSR835KFT; CD4 T cells) and GSE313371 (CD8 T cells). Thus, in the revised manuscript each of the three '
          'human datasets is compared with both mouse datasets, and the datasets of each species are compared with '
          'each other.'),
    ('b', 'Replaced the "splicing unchanged" category with formal equivalence testing (confidence intervals of '
          'ΔPSI and two one-sided tests), separating differential, equivalent and inconclusive splicing, and '
          'redefined high confidence events as those reproduced in every human-mouse dataset pair and absent from '
          'all within-species comparisons.'),
    ('b', 'Added effect-size analyses independent of significance thresholds: distributions of |ΔPSI| and '
          'correlations of PSI for human-mouse, human-human and mouse-mouse comparisons on the same events, and '
          'sensitivity analyses of the ΔPSI thresholds.'),
    ('b', 'Updated the exon orthology map to Ensembl release 116 and repeated all analyses with stricter exon '
          'identity thresholds (70% and 80%), with reciprocal best matches, and for events with identical human-mouse '
          'splice-site geometry.'),
    ('b', 'Compiled a sample-level metadata table and repeated the analysis with the most closely matched cell '
          'populations.'),
    ('b', 'Compared effect sizes between immune cells and fibroblasts on the same events and tested species x cell '
          'type interactions, and modeled the effect of expression and coverage on the classification and on the '
          'functional enrichment.'),
    ('p', 'Differences in age, sex and tissue between the human and mouse samples could explain the human-mouse '
          'differences. We will test for differences in splicing between male and female donors and between age '
          'groups in human datasets: male vs female donors of sorted naive B, CD4 and CD8 T cells, NK cells and '
          'monocytes (GSE107011, 2 vs 2 donors in each of these five cell types), and young (23-35 years) vs old '
          '(65-78 years) donors of naive CD4 T cells (GSE143213). In mouse, we will add sex (ENCODE spleen B cells, '
          '8 male vs 8 female mice) and age (naive CD4 T cells of 7-11-week vs 83-86-week-old mice, GSE184496, 6 vs '
          '6 mice) comparisons, and add them as supplementary data. The tissue aspect is more challenging, as most '
          'human datasets are from blood, and most mouse data are from tissues such as spleen, lymph nodes and bone '
          'marrow. We will compare the monocytes that were profiled in mouse blood (GSE116177) with the human blood '
          'monocytes, and look for additional samples from the same tissue that can be compared. As controls for '
          'tissue of origin within each species, we will compare NK cells from blood and spleen of the same human '
          'donors (GSE133383), NK cells from spleen and liver of the same mice (GSE233347), and mouse monocytes from '
          'bone marrow and blood (GSE116177).'),
    ('p', 'We therefore kindly request an extension of the revision deadline by four weeks, until [DATE]. We believe '
          'that these additional analyses will substantially strengthen the manuscript and fully address the '
          'reviewers\' concerns.'),
    ('p', 'Thank you for your consideration.'),
    ('p', 'Sincerely,'),
    ('p', 'Hadas Ner-Gaon, Naama Nudelman and Tal Shay'),
    ('p', 'Department of Life Sciences, Ben-Gurion University of the Negev, Beer-Sheva, Israel'),
]


def main():
    doc = Document()
    apply_style(doc)
    md = []
    for kind, text in PARAGRAPHS:
        if kind == 'h':
            doc.add_heading(text, level=1)
            md += [f'# {text}', '']
        elif kind == 'b':
            doc.add_paragraph(text, style='List Bullet')
            md += [f'- {text}']
        else:
            if md and md[-1].startswith('- '):
                md.append('')
            doc.add_paragraph(text)
            md += [text, '']
    finalize(doc)
    doc.save(OUT + '.docx')
    open(OUT + '.md', 'w').write('\n'.join(md) + '\n')
    print('wrote', OUT + '.docx')


if __name__ == '__main__':
    main()
