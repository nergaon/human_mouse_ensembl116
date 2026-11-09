"""House style for the Word files of the revision: Times New Roman 12, single line spacing (text and headings)."""
from docx.enum.text import WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, RGBColor


def apply_style(doc):
    for name in ('Normal', 'List Bullet', 'Title', 'Heading 1', 'Heading 2', 'Heading 3'):
        try:
            st = doc.styles[name]
        except KeyError:
            continue
        st.font.name = 'Times New Roman'
        st.element.get_or_add_rPr().get_or_add_rFonts().set(qn('w:eastAsia'), 'Times New Roman')
        st.font.size = Pt(12)
        if name.startswith('Heading') or name == 'Title':
            st.font.color.rgb = RGBColor(0, 0, 0)
        pf = st.paragraph_format
        pf.line_spacing_rule = WD_LINE_SPACING.SINGLE


def explicit_sizes(doc, half_points=24):
    """Write the size (sz and szCs) into the document defaults and into every run and paragraph mark that has none,
    so that Word shows 12 pt everywhere instead of falling back to the 11 pt default of the template."""
    def set_size(rpr):
        for tag in ('w:sz', 'w:szCs'):
            if rpr.find(qn(tag)) is None:
                el = OxmlElement(tag)
                el.set(qn('w:val'), str(half_points))
                rpr.append(el)
                # keep the schema order: sz/szCs before highlight, u, vertAlign, rtl, lang, ...
                later = [c for c in rpr if c.tag in AFTER_SZ]
                if later:
                    later[0].addprevious(el)
    defaults = doc.styles.element.find(qn('w:docDefaults'))
    rpr_default = defaults.find(qn('w:rPrDefault'))
    if rpr_default is not None:
        rpr = rpr_default.find(qn('w:rPr'))
        for tag in ('w:sz', 'w:szCs'):
            el = rpr.find(qn(tag))
            if el is not None:
                el.set(qn('w:val'), str(half_points))
        set_size(rpr)
    body = doc.element.body
    for r in body.iter(qn('w:r')):
        if r.find(qn('w:t')) is not None:
            set_size(r.get_or_add_rPr())
    for p in body.iter(qn('w:p')):
        set_size(_ppr_rpr(p))


def _ppr_rpr(p):
    ppr = p.get_or_add_pPr()
    rpr = ppr.find(qn('w:rPr'))
    if rpr is None:
        rpr = OxmlElement('w:rPr')
        # rPr is near the end of pPr: before sectPr and pPrChange only
        tail = [c for c in ppr if c.tag in (qn('w:sectPr'), qn('w:pPrChange'))]
        if tail:
            tail[0].addprevious(rpr)
        else:
            ppr.append(rpr)
    return rpr


AFTER_SZ = {qn('w:' + t) for t in ('highlight', 'u', 'effect', 'bdr', 'shd', 'fitText', 'vertAlign', 'rtl', 'cs',
                                    'em', 'lang', 'eastAsianLayout', 'specVanish', 'oMath')}


def explicit_fonts(doc, font='Times New Roman'):
    """Times New Roman for every script (ascii, hAnsi, eastAsia, cs) in the document defaults, in every style and in
    every run and paragraph mark; theme-font attributes are removed because they override the explicit fonts."""
    theme = [qn(f'w:{a}') for a in ('asciiTheme', 'hAnsiTheme', 'eastAsiaTheme', 'cstheme')]

    def set_font(rpr):
        rf = rpr.find(qn('w:rFonts'))
        if rf is None:
            rf = OxmlElement('w:rFonts')
            rstyle = rpr.find(qn('w:rStyle'))
            if rstyle is not None:
                rstyle.addnext(rf)
            else:
                rpr.insert(0, rf)
        for a in theme:
            rf.attrib.pop(a, None)
        for a in ('ascii', 'hAnsi', 'eastAsia', 'cs'):
            rf.set(qn(f'w:{a}'), font)

    styles = doc.styles.element
    rpr_default = styles.find(qn('w:docDefaults')).find(qn('w:rPrDefault'))
    set_font(rpr_default.find(qn('w:rPr')))
    for st in styles.iter(qn('w:style')):
        rpr = st.find(qn('w:rPr'))
        if rpr is None:
            if st.get(qn('w:type')) not in ('paragraph', 'character', 'table'):
                continue
            rpr = OxmlElement('w:rPr')
            tail = [c for c in st if c.tag in (qn('w:tblPr'), qn('w:trPr'), qn('w:tcPr'), qn('w:tblStylePr'))]
            if tail:
                tail[0].addprevious(rpr)
            else:
                st.append(rpr)
        set_font(rpr)
    for rf in styles.iter(qn('w:rFonts')):      # also the conditional formats of table styles
        set_font(rf.getparent())
    body = doc.element.body
    for r in body.iter(qn('w:r')):
        set_font(r.get_or_add_rPr())
    for p in body.iter(qn('w:p')):
        set_font(_ppr_rpr(p))


def finalize(doc):
    """Last step before saving: explicit Times New Roman and 12 pt everywhere, and a valid zoom element."""
    explicit_fonts(doc)
    explicit_sizes(doc)
    zoom = doc.settings.element.find(qn('w:zoom'))
    if zoom is not None and zoom.get(qn('w:percent')) is None:
        zoom.set(qn('w:percent'), '100')
