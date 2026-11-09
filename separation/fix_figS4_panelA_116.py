#!/usr/bin/env python3
"""
Figure S4 panel A (submitted Image_2.TIFF): the mouse junction coordinates pointed to mouse exons 4-6 of the
canonical transcript (orthologous to human exons 8-10), not to the orthologues of human exons 4-6. Replace them with the
orthologous mouse junctions used in the analysis (Table S7; Ensembl 116 transcript Puf60-202, ENSMUST00000100527,
whose exons 4, 5, 6 are 90, 51, 162 nt like human exons 4, 5, 6).

The digits are replaced glyph by glyph with digit glyphs taken from the same image (same font and size), centred on the
old glyph, so the layout is unchanged.
Output: ensembl116/revision_separation/figures/panels/PUF60_Image_2_fixed.tif
"""
import numpy as np
from PIL import Image

SRC = '/gpfs0/tals/users/nergaon/home/.claude/uploads/local_8a84ed4e-ab88-45b4-b7d9-0727a9d85737/Image_2.TIFF'
OUT = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_separation/figures/panels/PUF60_Image_2_fixed.tif'
X0 = 1700
# (top row, bottom row, old text, new text)
MOUSE = [(932, 967, '75,943,166–75,943,315', '75,944,417–75,946,456'),
         (986, 1021, '75,943,505–75,943,589', '75,946,506–75,947,391'),
         (1040, 1075, '75,943,166–75,943,589', '75,944,417–75,947,391')]
HUMAN = [(257, 292, '143,818,534–143,820,666'), (311, 346, '143,820,716–143,821,597'),
         (365, 400, '143,818,534–143,821,597')]


def glyphs(g, a, b):
    cols = np.where((g[a:b + 1, X0:] < 200).sum(0) > 0)[0] + X0
    segs, s, p = [], cols[0], cols[0]
    for c in cols[1:]:
        if c - p > 1:
            segs.append((s, p))
            s = c
        p = c
    segs.append((s, p))
    return segs


def main():
    im = np.asarray(Image.open(SRC).convert('RGB')).copy()
    g = np.asarray(Image.fromarray(im).convert('L'))
    lib = {}
    for a, b, text in HUMAN + [(a, b, old) for a, b, old, _ in MOUSE]:
        segs = glyphs(g, a, b)
        assert len(segs) >= len(text), (text, len(segs))
        for ch, (x1, x2) in zip(text, segs):
            if ch.isdigit():
                lib[ch] = im[a - 3:b + 4, x1:x2 + 1].copy()    # mouse lines last: their glyphs win
    for a, b, old, new in MOUSE:
        assert len(old) == len(new)
        segs = glyphs(g, a, b)
        for (x1, x2), co, cn in zip(segs, old, new):
            if co == cn:
                continue
            im[a - 3:b + 4, x1:x2 + 1] = 255
            gl = lib[cn]
            left = int(round((x1 + x2) / 2 - (gl.shape[1] - 1) / 2))
            region = im[a - 3:b + 4, left:left + gl.shape[1]]
            im[a - 3:b + 4, left:left + gl.shape[1]] = np.minimum(region, gl)
    Image.fromarray(im).save(OUT, compression='tiff_lzw', dpi=(300, 300))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
