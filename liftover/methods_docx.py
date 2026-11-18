#!/usr/bin/env python3
"""Methods text of the boundary-mapping correction (liftOver) and the read-based QC, as a Word file for review.
Output: ensembl116/liftover_check/Methods_boundary_mapping.docx"""
import sys

import docx

sys.path.insert(0, '/gpfs0/tals/projects/Analysis/human_mouse_exons/scripts/ensembl116/six')
from docx_style_116 import apply_style, explicit_sizes

OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/liftover_check/Methods_boundary_mapping.docx'

TEXT = [
    ('h', '2.1 Orthologous exon boundaries (added after the paragraph on the selection of the best orthologous '
          'boundary)'),
    ('p', 'Exon identity does not distinguish between annotated variants of the same exon whose boundaries are only a '
          'few nucleotides apart, for example an exon with two alternative 3\' splice sites. Therefore, for every '
          'orthologous boundary pair in which the human boundary or its mouse partner had another annotated splice site '
          'of the same type (donor or acceptor; all Ensembl 116 transcripts) within 30 nucleotides (65,919 of the '
          '390,686 pairs), the human boundary was projected onto the mouse genome with the UCSC hg38ToMm39 chain file '
          '(liftOver), and the projected position was compared with the mouse sites of the same type within 30 '
          'nucleotides of the selected partner. The projected position coincided with the selected partner in 63,080 '
          'pairs (95.7%). When it coincided, within two nucleotides, with another mouse site (1,050 pairs), the partner '
          'was replaced by that site only if it is a splice site of the Ensembl canonical transcript and the selected '
          'partner is not, and if it was not already the partner of another human boundary (459 pairs). The projection '
          'alone was not used to replace partners, because whole-genome alignments place gaps ambiguously near splice '
          'sites; the projection and the canonical transcript, both independent of the RNA-seq data, had to agree.'),
    ('h', '2.3 Defining human-mouse orthologous alternative splicing events (added at the end of the section)'),
    ('p', 'Splice-site quality control. An apparent difference between species can also arise when the two species '
          'use splice sites a few nucleotides apart and the event contains the site of one species and the unused, '
          'orthologous position in the other. To exclude such events, every junction of every event was compared with '
          'the junctions observed in the immune cell RNA-seq reads of each species. A junction was flagged if, in one '
          'species, it had fewer than 5% of the reads of a junction that is not part of the event, shares one end with '
          'it, has its other end within 30 nucleotides and has at least 100 reads, while in the other species it '
          'accounted for at least 10% of the reads of the event. Events with a flagged junction (21 of 6,810 events, 0.3%) were '
          'excluded from all analyses.'),
]


def main():
    d = docx.Document()
    apply_style(d)
    for kind, t in TEXT:
        p = d.add_paragraph()
        r = p.add_run(t)
        r.bold = kind == 'h'
    explicit_sizes(d)
    d.save(OUT)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
