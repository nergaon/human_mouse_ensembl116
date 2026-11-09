#!/usr/bin/env python3
"""
Gene figures of the revision: panel A (exon diagram) from the submitted figure, panel B regenerated
from the Ensembl 116 data (gene_panel_116.py), on an A4 page at 300 dpi. New figures (DRAWN, and the candidates)
have panel A drawn from the Ensembl 116 data too (gene_diagram_116.py).

Panel A is everything above the 'B' panel label of the submitted figure (found in the left margin).
Output: ensembl116/revision_separation/figures/<Figure>.tif / .pdf
"""
import os
import subprocess

import numpy as np
from PIL import Image, ImageDraw, ImageFont

UP = '/gpfs0/tals/users/nergaon/home/.claude/uploads/local_8a84ed4e-ab88-45b4-b7d9-0727a9d85737/'
FIG = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/revision_separation/figures/'
PANELS = FIG + 'panels/'
FONT = '/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf'
FONT_REG = '/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf'
A4 = (2480, 3508)

# new figure name: (source of panel A, new panel B image)
FIGURES = {
    'Figure_4': (UP + 'figure_3_v11_USP16.pdf', PANELS + 'USP16_clu_52144.png'),
    'Figure_S4': (PANELS + 'PUF60_Image_2_fixed.tif', PANELS + 'PUF60_clu_85337.png'),  # mouse coordinates fixed (fix_figS4_panelA_116.py)
    'Figure_S5': (UP + 'Image_3.TIFF', PANELS + 'RPL22L1_clu_64510.png'),
    'Figure_S8': (UP + 'DOCK8_figure_5_v6.pdf', PANELS + 'DOCK8_clu_90159.png'),
    'Figure_S9': (UP + 'Image_6.TIFF', PANELS + 'PTPRC_clu_44342.png'),
    'Figure_S10': (UP + 'Image_7.TIFF', PANELS + 'CD47_clu_62548.png'),
    'Figure_S6': (UP + 'Image_8.TIFF', PANELS + 'UCKL1_clu_51932.png'),
}


def load_page(path):
    if path.lower().endswith('.pdf'):
        out = FIG + 'tmp_page'
        subprocess.run(['pdftoppm', '-r', '300', '-png', '-f', '1', '-l', '1', '-singlefile', path, out], check=True)
        im = Image.open(out + '.png').convert('RGB')
        os.remove(out + '.png')
    else:
        im = Image.open(path).convert('RGB')
    return im.resize(A4) if abs(im.size[0] - A4[0]) > 5 else im


def panel_b_row(im):
    """Top row of the 'B' label: second group of dark pixels in the left margin."""
    g = np.asarray(im.convert('L'))
    h, w = g.shape
    rows = np.where((g[:, int(0.01 * w):int(0.06 * w)] < 128).sum(1) > 0)[0]
    groups, start, prev = [], None, None
    for r in rows:
        if start is None:
            start = prev = r
        elif r - prev <= 5:
            prev = r
        else:
            groups.append((start, prev))
            start = prev = r
    groups.append((start, prev))
    return groups[1][0]


# figures whose left margin holds other labels (row names): 'B' label row set by hand (300 dpi)
B_ROW = {'Figure_4': 1469}


def assemble(name, src, panel_b):
    page = load_page(src)
    cut = B_ROW.get(name, panel_b_row(page)) - 15
    panel_a = page.crop((0, 0, page.size[0], cut))
    b = Image.open(panel_b).convert('RGB')
    width = int(A4[0] * 0.92)
    b = b.resize((width, int(b.size[1] * width / b.size[0])), Image.LANCZOS)
    out = Image.new('RGB', A4, 'white')
    out.paste(panel_a, (0, 0))
    y_b = cut + 40
    out.paste(b, ((A4[0] - width) // 2 + 30, y_b + 50))
    draw = ImageDraw.Draw(out)
    draw.text((int(0.013 * A4[0]), y_b), 'B', fill='black', font=ImageFont.truetype(FONT, 44))
    label = name.replace('_', ' ')
    draw.text((int(0.03 * A4[0]), min(y_b + 50 + b.size[1] + 120, A4[1] - 120)), label, fill='black',
              font=ImageFont.truetype(FONT_REG, 44))
    out.save(FIG + name + '.tif', compression='tiff_lzw', dpi=(300, 300))
    out.save(FIG + name + '.pdf', resolution=300)
    print('wrote', name, 'panel A height', cut)


# figures with panel A drawn by gene_diagram_116.py (Figure 5 was HNRNPK in the submitted manuscript)
DRAWN = {
    'Figure_5': (PANELS + 'SRSF7_clu_60271_A.png', PANELS + 'SRSF7_clu_60271.png'),
}
# candidate figures, not in the manuscript (FERMT3, removed 2026-11-02, was a boundary-mapping artifact)
CANDIDATES = {
    # DS in nearly all H-M and also all H-H pairs (dataset-driven DS example; user asked 2026-11-03)
    'IKZF1_clu_83293': (PANELS + 'IKZF1_clu_83293_A.png', PANELS + 'IKZF1_clu_83293.png'),
    # top species-separated event (species_separation_116.py, gap 0.83 in all five cell types)
    'TRAF3IP3_clu_44912': (PANELS + 'TRAF3IP3_clu_44912_A.png', PANELS + 'TRAF3IP3_clu_44912.png'),
}


def assemble_drawn(name, panel_a, panel_b, folder=FIG + 'candidates/', label=False):
    out = Image.new('RGB', A4, 'white')
    draw = ImageDraw.Draw(out)
    width = int(A4[0] * 0.92)
    y = 80
    for letter, path in (('A', panel_a), ('B', panel_b)):
        im = Image.open(path).convert('RGB')
        w = int(width * (0.8 if letter == 'A' else 1))
        im = im.resize((w, int(im.size[1] * w / im.size[0])), Image.LANCZOS)
        draw.text((int(0.013 * A4[0]), y), letter, fill='black', font=ImageFont.truetype(FONT, 44))
        out.paste(im, ((A4[0] - w) // 2 + 30, y + 50))
        y += 50 + im.size[1] + 80
    if label:
        draw.text((int(0.03 * A4[0]), min(y + 40, A4[1] - 120)), name.replace('_', ' '), fill='black',
                  font=ImageFont.truetype(FONT_REG, 44))
    os.makedirs(folder, exist_ok=True)
    out.save(folder + f'{name}.tif', compression='tiff_lzw', dpi=(300, 300))
    out.save(folder + f'{name}.pdf', resolution=300)
    print('wrote', name)


def main():
    for name, (a, b) in CANDIDATES.items():
        assemble_drawn(name, a, b)
    for name, (a, b) in DRAWN.items():
        assemble_drawn(name, a, b, folder=FIG, label=True)
    for name, (src, b) in FIGURES.items():
        if not os.path.exists(src):
            print('missing panel A source for', name, src)
            continue
        assemble(name, src, b)


if __name__ == '__main__':
    main()
