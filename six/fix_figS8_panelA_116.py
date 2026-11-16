#!/usr/bin/env python3
"""
Figure S8 panel A (submitted DOCK8_figure_5_v6.pdf): the 'H' of 'Human DOCK8' is hidden behind the panel label 'A'.
Redraw the title 'Human DOCK8' left-aligned with 'Mouse Dock8' (same font, Liberation Sans Bold, size fitted to the
'Mouse Dock8' text of the same page).
Output: ensembl116/revision_6ds/figures/panels/DOCK8_panelA_fixed.png (page at 300 dpi)
"""
import numpy as np
from PIL import Image, ImageDraw, ImageFont

import assemble_gene_figures_116 as a

OUT = a.PANELS + 'DOCK8_panelA_fixed.png'


def bbox(g, y0, y1, x0, x1):
    reg = g[y0:y1, x0:x1] < 128
    rows, cols = np.where(reg.any(1))[0], np.where(reg.any(0))[0]
    return x0 + cols.min(), y0 + rows.min(), x0 + cols.max(), y0 + rows.max()


def main():
    page = a.load_page(a.UP + 'DOCK8_figure_5_v6.pdf')
    g = np.asarray(page.convert('L'))
    mx0, my0, mx1, my1 = bbox(g, 560, 660, 0, 700)            # 'Mouse Dock8'
    hx0, hy0, hx1, hy1 = bbox(g, 30, 130, 200, 700)           # 'Human DOCK8' (visible part, right of the 'A')
    # font size whose rendering of 'Mouse Dock8' has the measured width
    size = min(range(20, 80), key=lambda s: abs(ImageFont.truetype(a.FONT, s).getbbox('Mouse Dock8')[2] -
                                                ImageFont.truetype(a.FONT, s).getbbox('Mouse Dock8')[0] - (mx1 - mx0)))
    font = ImageFont.truetype(a.FONT, size)
    d = ImageDraw.Draw(page)
    d.rectangle((hx0 - 2, hy0 - 4, hx1 + 4, hy1 + 4), fill='white')
    l, t, _, b = font.getbbox('Human DOCK8')
    d.text((mx0 - l, hy1 - b), 'Human DOCK8', fill='black', font=font)   # same left edge as Mouse, same baseline
    page.save(OUT)
    print('wrote', OUT, 'font size', size)


if __name__ == '__main__':
    main()
