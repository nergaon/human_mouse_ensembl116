#!/usr/bin/env python3
"""
Build the revised manuscript as tracked changes against the submitted HN44.docx, plus a clean copy.

- changed paragraphs: word-level diff of old vs new text; unchanged words keep their original runs
  (formatting), deleted words become <w:del>, new words <w:ins> with the formatting of the preceding run
- citations: all numbered citations, old and new ({key}), are renumbered by order of first appearance
  (title -> figure legends); the reference list is rebuilt in that order
- Table 1: the old table is marked deleted and the new table (table1.tsv, facts_116.py) inserted
Usage: make_revised_docx_116.py <HN44.docx> <unpacked + run-merged HN44 dir> <output dir>
"""
import datetime
import html
import os
import re
import shutil
import subprocess
import sys
from difflib import SequenceMatcher

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from revised_text_116 import EDITS, INSERTS, NEW_REFERENCES  # noqa: E402

AUTHOR = 'Hadas Ner-Gaon'
DATE = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%SZ')
REF_FIRST, REF_LAST = 234, 295          # reference list paragraphs in HN44
TABLE_FIRST, TABLE_LAST = 105, 219      # paragraphs inside Table 1
ED = '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/'
_id = [9000]


def nid():
    _id[0] += 1
    return _id[0]


def esc(t):
    return html.escape(t, quote=False)


# ------------------------------------------------------------------ paragraph parsing
P_RE = re.compile(r'<w:p[ >].*?</w:p>', re.S)


def runs_of(p):
    """List of (char, rPr) for the text of a paragraph; tabs as '\\t'."""
    chars = []
    for r in re.finditer(r'<w:r[ >].*?</w:r>', p, re.S):
        rx = r.group(0)
        rpr = re.search(r'<w:rPr>.*?</w:rPr>', rx, re.S)
        rpr = rpr.group(0) if rpr else ''
        for item in re.finditer(r'<w:t(?: [^>]*)?>([^<]*)</w:t>|<w:tab/>', rx):
            if item.group(0) == '<w:tab/>':
                chars.append(('\t', rpr))
            else:
                chars.extend((c, rpr) for c in html.unescape(item.group(1)))
    return chars


def text_of(p):
    return ''.join(c for c, _ in runs_of(p))


def run_xml(text, rpr, kind='normal'):
    out = []
    for part in re.split(r'(\t)', text):
        if not part:
            continue
        if part == '\t':
            out.append('<w:tab/>')
        elif kind == 'del':
            out.append(f'<w:delText xml:space="preserve">{esc(part)}</w:delText>')
        else:
            out.append(f'<w:t xml:space="preserve">{esc(part)}</w:t>')
    body = f'<w:r>{rpr}{"".join(out)}</w:r>'
    if kind == 'del':
        return f'<w:del w:id="{nid()}" w:author="{AUTHOR}" w:date="{DATE}">{body}</w:del>'
    if kind == 'ins':
        return f'<w:ins w:id="{nid()}" w:author="{AUTHOR}" w:date="{DATE}">{body}</w:ins>'
    return body


def grouped(chars, kind):
    out, buf, cur = [], '', None
    for c, r in chars:
        if r != cur and buf:
            out.append(run_xml(buf, cur, kind))
            buf = ''
        cur = r
        buf += c
    if buf:
        out.append(run_xml(buf, cur, kind))
    return ''.join(out)


TOKEN = re.compile(r'\s+|\w+|[^\w\s]', re.U)


def redline_paragraph(p, new_text):
    chars = runs_of(p)
    old = ''.join(c for c, _ in chars)
    if old == new_text:
        return p
    a = [m.group(0) for m in TOKEN.finditer(old)]
    b = [m.group(0) for m in TOKEN.finditer(new_text)]
    pos_a = [0]
    for t in a:
        pos_a.append(pos_a[-1] + len(t))
    body = []
    default_rpr = chars[0][1] if chars else ''
    for op, i1, i2, j1, j2 in SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        seg = chars[pos_a[i1]:pos_a[i2]]
        if op == 'equal':
            body.append(grouped(seg, 'normal'))
            continue
        if op in ('delete', 'replace'):
            body.append(grouped(seg, 'del'))
        if op in ('insert', 'replace'):
            rpr = chars[pos_a[i1] - 1][1] if pos_a[i1] > 0 else default_rpr
            body.append(run_xml(''.join(b[j1:j2]), rpr, 'ins'))
    ppr = re.search(r'<w:pPr>.*?</w:pPr>', p, re.S)
    head = re.match(r'<w:p[^>]*>', p).group(0)
    return head + (ppr.group(0) if ppr else '') + ''.join(body) + '</w:p>'


def inserted_paragraph(template, text):
    ppr = re.search(r'<w:pPr>.*?</w:pPr>', template, re.S)
    ppr = ppr.group(0) if ppr else '<w:pPr></w:pPr>'
    mark = f'<w:ins w:id="{nid()}" w:author="{AUTHOR}" w:date="{DATE}"/>'
    if '<w:rPr>' in ppr:
        ppr = ppr.replace('<w:rPr>', '<w:rPr>' + mark, 1)
    else:
        ppr = ppr.replace('</w:pPr>', f'<w:rPr>{mark}</w:rPr></w:pPr>')
    chars = runs_of(template)
    rpr = chars[0][1] if chars else ''
    return f'<w:p>{ppr}{run_xml(text, rpr, "ins") if text else ""}</w:p>'


def deleted_paragraph(p):
    chars = runs_of(p)
    ppr = re.search(r'<w:pPr>.*?</w:pPr>', p, re.S)
    ppr = ppr.group(0) if ppr else '<w:pPr></w:pPr>'
    mark = f'<w:del w:id="{nid()}" w:author="{AUTHOR}" w:date="{DATE}"/>'
    if '<w:rPr>' in ppr:
        ppr = ppr.replace('<w:rPr>', '<w:rPr>' + mark, 1)
    else:
        ppr = ppr.replace('</w:pPr>', f'<w:rPr>{mark}</w:rPr></w:pPr>')
    head = re.match(r'<w:p[^>]*>', p).group(0)
    return head + ppr + grouped(chars, 'del') + '</w:p>'


# ------------------------------------------------------------------ citations
CIT = re.compile(r'\(((?:\d+|\{[a-z]+\})(?:\s*[,\-–]\s*(?:\d+|\{[a-z]+\}))*)\)')


def parse_group(g):
    items = []
    for part in re.split(r'\s*,\s*', g):
        m = re.fullmatch(r'(\d+)\s*[\-–]\s*(\d+)', part)
        if m:
            items.extend(range(int(m.group(1)), int(m.group(2)) + 1))
        elif part.isdigit():
            items.append(int(part))
        else:
            items.append(part.strip('{}'))
    return items


def compress(nums):
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(f'{nums[i]}-{nums[j]}' if j - i >= 2 else ', '.join(map(str, nums[i:j + 1])))
        i = j + 1
    return ', '.join(out)


def renumber(texts, n_old):
    order = {}
    for t in texts:
        for m in CIT.finditer(t):
            items = parse_group(m.group(1))
            if any(isinstance(x, int) and x > n_old for x in items):
                continue
            for x in items:
                order.setdefault(x, len(order) + 1)

    def sub(m):
        items = parse_group(m.group(1))
        if any(isinstance(x, int) and x > n_old for x in items):
            return m.group(0)
        return '(' + compress([order[x] for x in items]) + ')'
    return [CIT.sub(sub, t) for t in texts], order


# ------------------------------------------------------------------ Table 1
CELLS = ['CD4T', 'CD8T', 'NveB', 'NK', 'Mono']


def table1_rows():
    t = pd.read_csv(ED + 'revision_files/table1.tsv', sep='\t').set_index(['kind', 'cell'])
    hc = pd.read_csv(ED + 'events/id50/classification/summary_by_cell.tsv', sep='\t', index_col=0)
    st = pd.read_csv(ED + 'events/id50/classification/event_status_by_cell.tsv', sep='\t')
    tc = st[st['status_high'] != 'not testable'].groupby('cell').size()
    fib = t.loc[('fibro', 'Fibro')]

    def row(k, col, fibro=''):
        return [t.loc[(k, c), col] if (k, c) in t.index else '' for c in CELLS] + [fibro]

    def hcrow(col):
        return [f'{hc.loc[c, col]} ({100 * hc.loc[c, col] / tc[c]:.1f}%)' for c in CELLS] + ['']
    R = [('header', ['Comparison category', 'CD4T', 'CD8T', 'B', 'NK', 'Mono', 'Fibroblasts'])]
    R.append(('section', ['Human-Mouse alternative splicing events (immune cells: 6 dataset pairs, classified per '
                          'event from two pairs that share no dataset; fibroblasts: 1 pair)']))
    # fibroblasts: one human-mouse pair, so the event status of that pair, without confidence levels
    fs = pd.read_csv(ED + 'events/id50/effects/event_effects.tsv', sep='\t')
    fs = fs.loc[(fs['comparison'] == 'FH_FM') & (fs['status_m0.1'] != 'not_testable'), 'status_m0.1'].value_counts()
    fn = int(fs.sum())

    def fcell(k):
        return f'{fs.get(k, 0)} ({100 * fs.get(k, 0) / fn:.1f}%)'
    R.append(('row', ['Testable events'] + [str(tc[c]) for c in CELLS] + [str(fn)]))
    R.append(('row', ['Differentially spliced events, high confidence'] + hcrow('high_divergent')[:-1] + ['–']))
    R.append(('row', ['Differentially spliced events, high + low confidence'] + hcrow('high_low_divergent')[:-1]
              + [fcell('DS')]))
    R.append(('row', ['Equivalent events, high confidence'] + hcrow('high_conserved')[:-1] + ['–']))
    R.append(('row', ['Equivalent events, high + low confidence'] + hcrow('high_low_conserved')[:-1]
              + [fcell('equivalent')]))
    R.append(('row', ['Inconclusive events'] + hcrow('high_low_inconclusive')[:-1] + [fcell('inconclusive')]))
    R.append(('section', ['Human-Mouse, matched populations (M1 naive T cells, blood monocytes; 3 pairs)']))
    R.append(('row', ['Differentially spliced events, %'] + row('matched_cross', 'DS')))
    R.append(('section', ['Human-Human alternative splicing events (3 dataset pairs)']))
    R.append(('row', ['Comparable events'] + row('HH', 'n')))
    R.append(('row', ['Differentially spliced events [p≤0.05 & max(|ΔPSI|)≥0.1], %'] + row('HH', 'DS')))
    R.append(('row', ['Equivalent events [90% CI of ΔPSI within ±0.1], %'] + row('HH', 'EQ')))
    R.append(('row', ['Inconclusive events, %'] + row('HH', 'INC')))
    R.append(('section', ['Mouse-Mouse alternative splicing events (1 dataset pair)']))
    R.append(('row', ['Comparable events'] + row('MM', 'n')))
    R.append(('row', ['Differentially spliced events [p≤0.05 & max(|ΔPSI|)≥0.1], %'] + row('MM', 'DS')))
    R.append(('row', ['Equivalent events [90% CI of ΔPSI within ±0.1], %'] + row('MM', 'EQ')))
    R.append(('row', ['Inconclusive events, %'] + row('MM', 'INC')))
    return R


def table1_xml():
    widths = [3402] + [945] * 6
    border = ('<w:tcBorders><w:top w:val="single" w:sz="6" w:space="0" w:color="auto"/><w:left w:val="single" '
              'w:sz="6" w:space="0" w:color="auto"/><w:bottom w:val="single" w:sz="6" w:space="0" w:color="auto"/>'
              '<w:right w:val="single" w:sz="6" w:space="0" w:color="auto"/></w:tcBorders>')
    mar = ('<w:tcMar><w:top w:w="40" w:type="dxa"/><w:left w:w="80" w:type="dxa"/><w:bottom w:w="40" '
           'w:type="dxa"/><w:right w:w="80" w:type="dxa"/></w:tcMar>')

    def cell(text, w, span=1, bold=False):
        rpr = ('<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="Times New Roman" '
               'w:cs="Times New Roman"/>' + ('<w:b/><w:bCs/>' if bold else '')
               + '<w:sz w:val="24"/><w:szCs w:val="24"/></w:rPr>')
        gs = f'<w:gridSpan w:val="{span}"/>' if span > 1 else ''
        mark = f'<w:ins w:id="{nid()}" w:author="{AUTHOR}" w:date="{DATE}"/>'
        # explicit left-to-right: the template's Normal style is right-to-left (bidi)
        p = (f'<w:p><w:pPr><w:bidi w:val="0"/><w:spacing w:after="0" w:line="240" w:lineRule="auto"/>'
             f'<w:jc w:val="left"/><w:rPr>{mark}{TNR_FONTS}<w:sz w:val="24"/><w:szCs w:val="24"/></w:rPr></w:pPr>'
             + (run_xml(text, rpr, 'ins') if text else '') + '</w:p>')
        return f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/>{gs}{border}{mar}<w:vAlign w:val="center"/></w:tcPr>{p}</w:tc>'
    rows = []
    for kind_, vals in table1_rows():
        trpr = f'<w:trPr><w:cantSplit/><w:ins w:id="{nid()}" w:author="{AUTHOR}" w:date="{DATE}"/></w:trPr>'
        if kind_ == 'section':
            cells = cell(vals[0], sum(widths), span=7, bold=True)
        else:
            cells = ''.join(cell(str(v), w, bold=(kind_ == 'header')) for v, w in zip(vals, widths))
        rows.append(f'<w:tr>{trpr}{cells}</w:tr>')
    grid = ''.join(f'<w:gridCol w:w="{w}"/>' for w in widths)
    return (f'<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:bidiVisual w:val="0"/><w:tblW w:w="{sum(widths)}" w:type="dxa"/><w:tblLayout w:type="fixed"/>'
            f'<w:tblLook w:val="04A0" w:firstRow="1" w:lastRow="0" w:firstColumn="1" w:lastColumn="0" '
            f'w:noHBand="0" w:noVBand="1"/></w:tblPr><w:tblGrid>{grid}</w:tblGrid>{"".join(rows)}</w:tbl>')


def delete_table(tbl):
    tbl = P_RE.sub(lambda m: deleted_paragraph(m.group(0)), tbl)
    return re.sub(r'</w:trPr>', lambda m: f'<w:del w:id="{nid()}" w:author="{AUTHOR}" w:date="{DATE}"/></w:trPr>',
                  re.sub(r'<w:tr( [^>]*)?>(?!<w:trPr>)', lambda m: m.group(0) + '<w:trPr></w:trPr>', tbl))


# ------------------------------------------------------------------ clean copy
def accept_all(xml):
    """Accept the tracked changes made by this script (HN44 itself had none)."""
    xml = re.sub(r'<w:del w:id="\d+" w:author="[^"]*" w:date="[^"]*">.*?</w:del>', '', xml, flags=re.S)
    # paragraphs whose mark is deleted were deleted entirely by this script
    xml = re.sub(r'<w:p[ >](?:(?!</w:p>).)*?<w:rPr><w:del w:id="\d+"[^>]*/>.*?</w:p>', '', xml, flags=re.S)
    xml = re.sub(r'<w:tr[ >](?:(?!</w:tr>).)*?<w:del w:id="\d+"[^/]*/></w:trPr>.*?</w:tr>', '', xml, flags=re.S)
    xml = re.sub(r'<w:tbl>(?:(?!<w:tr[ >]).)*?</w:tbl>', '', xml, flags=re.S)  # tables left without rows
    xml = re.sub(r'<w:ins w:id="\d+" w:author="[^"]*" w:date="[^"]*"/>', '', xml)
    xml = re.sub(r'<w:ins w:id="\d+" w:author="[^"]*" w:date="[^"]*">(.*?)</w:ins>', r'\1', xml, flags=re.S)
    return xml


TNR_FONTS = ('<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="Times New Roman" '
             'w:cs="Times New Roman"/>')
SINGLE = '<w:spacing w:before="0" w:beforeAutospacing="0" w:after="0" w:afterAutospacing="0" w:line="240" w:lineRule="auto"/>'


def house_style(work):
    """Times New Roman 12 pt, single line spacing and no space before / after paragraphs (title stays 16 pt;
    headings keep their space before)."""
    tnr = 'Times New Roman'
    path = os.path.join(work, 'word', 'styles.xml')
    st = open(path, encoding='utf8').read()
    st = re.sub(r'(w:(?:ascii|hAnsi|cs|eastAsia))="(?:Calibri|Calibri Light|Tahoma|Helvetica|Arial)"', rf'\1="{tnr}"', st)
    st = re.sub(r'w:line="\d+" w:lineRule="auto"', 'w:line="240" w:lineRule="auto"', st)

    def flat(m):  # paragraph spacing of every non-heading style (character spacing has w:val and is kept)
        style = m.group(0)
        if re.search(r'w:styleId="(Heading\d|Title)"', style):
            return style
        return re.sub(r'<w:spacing (?![^>]*w:val=)[^>]*/>', '<w:snapToGrid w:val="0"/>' + SINGLE, style)
    st = re.sub(r'<w:style [^>]*>.*?</w:style>', flat, st, flags=re.S)
    st = re.sub(r'(<w:pPrDefault><w:pPr>)<w:spacing [^>]*/>', rf'\1{SINGLE}', st)
    st = re.sub(r'<w:rFonts [^>]*/>', TNR_FONTS, st)  # explicit fonts everywhere (no theme or fallback fonts)
    st = re.sub(r'(<w:docDefaults>.*?)<w:sz w:val="\d+"/><w:szCs w:val="\d+"/>',
                r'\1<w:sz w:val="24"/><w:szCs w:val="24"/>', st, count=1, flags=re.S)
    open(path, 'w', encoding='utf8').write(st)
    path = os.path.join(work, 'word', 'theme', 'theme1.xml')
    th = open(path, encoding='utf8').read()
    th = re.sub(r'(<a:(?:major|minor)Font>\s*<a:latin typeface=")[^"]+"', rf'\1{tnr}"', th)
    th = re.sub(r'<a:(ea|cs) typeface="[^"]*"', rf'<a:\1 typeface="{tnr}"', th)
    open(path, 'w', encoding='utf8').write(th)


def restyle_document(xml):
    xml = re.sub(r'w:ascii="Helvetica" w:hAnsi="Helvetica"', 'w:ascii="Times New Roman" w:hAnsi="Times New Roman"', xml)
    xml = re.sub(r'<w:rFonts [^>]*/>', TNR_FONTS, xml)  # same font metrics in every run
    xml = xml.replace('<w:sz w:val="22"/>', '<w:sz w:val="24"/>').replace('<w:szCs w:val="22"/>', '<w:szCs w:val="24"/>')
    return uniform_paragraphs(explicit_sizes(xml))


RPR_AFTER_SZ = re.compile(r'<w:(?:highlight|u|effect|bdr|shd|fitText|vertAlign|rtl|cs|em|lang|eastAsianLayout|'
                          r'specVanish|oMath)[ />]')


def explicit_sizes(xml):
    """12 pt written in every run (and paragraph mark) that has no size of its own, so that no viewer falls back
    to another size; runs with their own size and comment/note references are left alone."""
    size = '<w:sz w:val="24"/><w:szCs w:val="24"/>'

    def fix_rpr(m):
        inner = m.group(1)
        if '<w:sz ' in inner or re.search(r'<w:rStyle w:val="\w*(?:Comment|Footnote|Endnote)', inner):
            return m.group(0)
        k = RPR_AFTER_SZ.search(inner)
        pos = k.start() if k else len(inner)
        return '<w:rPr>' + inner[:pos] + size + inner[pos:] + '</w:rPr>'

    xml = re.sub(r'<w:rPr>((?:(?!</w:rPr>).)*?)</w:rPr>', fix_rpr, xml, flags=re.S)
    xml = re.sub(r'<w:rPr/>', f'<w:rPr>{size}</w:rPr>', xml)
    return re.sub(r'(<w:r(?: [^>]*)?>)(?!<w:rPr>)', rf'\1<w:rPr>{size}</w:rPr>', xml)


PPR_AFTER = ('ind', 'contextualSpacing', 'mirrorIndents', 'suppressOverlap', 'jc', 'textDirection', 'textAlignment',
             'textboxTightWrap', 'outlineLvl', 'divId', 'cnfStyle', 'rPr', 'sectPr', 'pPrChange')


def uniform_paragraphs(xml):
    """Every paragraph gets the same explicit line layout: single spacing, 0 pt before/after, no snapping to the
    document grid (the template defines an 18 pt grid), no contextual spacing."""
    layout = '<w:snapToGrid w:val="0"/>' + SINGLE

    def fix_ppr(m):
        inner = m.group(1)
        inner = re.sub(r'<w:(?:snapToGrid|contextualSpacing)(?: [^>]*)?/>', '', inner)
        inner = re.sub(r'<w:spacing (?![^>]*w:val=)[^>]*/>', '', inner)
        pos = min([k for k in (inner.find(f'<w:{tag}') for tag in PPR_AFTER) if k != -1] + [len(inner)])
        own = layout
        if re.search(r'<w:pStyle w:val="(Heading\d|Title)"', inner):  # headings keep their space before
            own = '<w:snapToGrid w:val="0"/><w:spacing w:line="240" w:lineRule="auto"/>'
        return '<w:pPr>' + inner[:pos] + own + inner[pos:] + '</w:pPr>'

    xml = re.sub(r'<w:pPr>((?:(?!</w:pPr>).)*?)</w:pPr>', fix_ppr, xml, flags=re.S)
    # paragraphs without properties
    xml = re.sub(r'(<w:p(?: [^>]*)?>)(?!<w:pPr>)', rf'\1<w:pPr>{layout}</w:pPr>', xml)
    return xml



# ------------------------------------------------------------------ main
def main():
    src_docx, unpacked, out_dir = sys.argv[1:4]
    os.makedirs(out_dir, exist_ok=True)
    work = os.path.join(out_dir, 'work')
    if os.path.exists(work):
        shutil.rmtree(work)
    shutil.copytree(unpacked, work)
    doc_path = os.path.join(work, 'word', 'document.xml')
    xml = open(doc_path, encoding='utf8').read()
    body_start = xml.find('<w:body>')
    spans = [(m.start() + body_start, m.end() + body_start) for m in P_RE.finditer(xml[body_start:])]
    paras = [xml[s:e] for s, e in spans]
    old_text = [text_of(p) for p in paras]

    # old references
    old_refs = {}
    for i in range(REF_FIRST, REF_LAST + 1):
        m = re.match(r'\((\d+)\)\s*(.*)', old_text[i])
        old_refs[int(m.group(1))] = m.group(2)
    n_old = max(old_refs)

    # new text of every body paragraph (references and table excluded), in document order, incl. inserts
    seq = []  # (kind, index, text): kind 'para' or 'ins' (index = paragraph after which it is inserted)
    for i in range(len(paras)):
        if REF_FIRST <= i <= REF_LAST or TABLE_FIRST <= i <= TABLE_LAST:
            continue
        seq.append(('para', i, EDITS.get(i, old_text[i])))
        for k, (tmpl, text) in enumerate(INSERTS.get(i, [])):
            seq.append(('ins', (i, k, tmpl), text))
    texts, order = renumber([t for _, _, t in seq], n_old)
    seq = [(k, i, t) for (k, i, _), t in zip(seq, texts)]
    new_refs = sorted(order.items(), key=lambda kv: kv[1])
    dropped = sorted(set(old_refs) - {k for k in order if isinstance(k, int)})
    print('references:', len(new_refs), 'dropped old refs:', dropped)

    # build replacements (from the end, so spans stay valid)
    repl = {}
    inserts_after = {}
    for kind_, idx, text in seq:
        if kind_ == 'para':
            if text != old_text[idx]:
                repl[idx] = redline_paragraph(paras[idx], text)
        else:
            i, k, tmpl = idx
            inserts_after.setdefault(i, []).append(inserted_paragraph(paras[tmpl], text))
    # reference list: delete old, insert new after the 'References' heading (REF_FIRST - 1)
    for i in range(REF_FIRST, REF_LAST + 1):
        repl[i] = deleted_paragraph(paras[i])
    ref_pars = []
    for key, num in new_refs:
        text = old_refs[key] if isinstance(key, int) else NEW_REFERENCES[key]
        ref_pars.append(inserted_paragraph(paras[REF_FIRST], f'({num}) {text}'))
    inserts_after.setdefault(REF_FIRST - 1, []).extend(ref_pars)

    out = []
    last = 0
    for i, (s, e) in enumerate(spans):
        if TABLE_FIRST <= i <= TABLE_LAST:
            continue  # handled with the table below
        out.append(xml[last:s])
        out.append(repl.get(i, paras[i]))
        out.extend(inserts_after.get(i, []))
        last = e
    out.append(xml[last:])
    new_xml = ''.join(out)
    # Table 1: the old table text between the caption and the next paragraph was skipped above; rebuild it
    tbl_old = re.search(r'<w:tbl>.*?</w:tbl>', xml, re.S).group(0)
    new_xml = new_xml.replace(tbl_old, delete_table(tbl_old) + '<w:p/>' + table1_xml(), 1) \
        if tbl_old in new_xml else new_xml
    new_xml = restyle_document(new_xml)
    open(doc_path, 'w', encoding='utf8').write(new_xml)
    house_style(work)

    tracked = os.path.join(out_dir, 'Human_Mouse_Splicing_revised_tracked.docx')
    if os.path.exists(tracked):
        os.remove(tracked)
    subprocess.run(['zip', '-q', '-X', '-r', os.path.abspath(tracked), '.'], cwd=work, check=True)
    print('wrote', tracked, '| changed paragraphs:', len(repl) - (REF_LAST - REF_FIRST + 1),
          '| inserted paragraphs:', sum(len(v) for v in inserts_after.values()))
    open(doc_path, 'w', encoding='utf8').write(accept_all(new_xml))
    clean = os.path.join(out_dir, 'Human_Mouse_Splicing_revised_clean.docx')
    if os.path.exists(clean):
        os.remove(clean)
    subprocess.run(['zip', '-q', '-X', '-r', os.path.abspath(clean), '.'], cwd=work, check=True)
    print('wrote', clean)
    with open(os.path.join(out_dir, 'reference_map.tsv'), 'w') as f:
        f.write('new\told_or_key\n')
        for key, num in new_refs:
            f.write(f'{num}\t{key}\n')


if __name__ == '__main__':
    main()
