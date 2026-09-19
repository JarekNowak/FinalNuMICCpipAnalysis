"""prl_wordcount.py -- word count of a revtex Letter against the PRL limit, APS rules.

PRL counts the body text (not abstract, authors, affiliations or references) plus word-equivalents
for figures and tables: a single-column figure is 150/aspect + 20 words, a double-column one
300/(0.5*aspect) + 40, where aspect = width/height; a single-column table is 13 + 6.5 words per
line and a double-column one 26 + 13 per line. Equations are counted as text. The limit is 3750.

    python3 report/tools/prl_wordcount.py report/paper_prl_proton_tagged.tex
"""
import re, sys, os, subprocess

LIMIT = 3750


def aspect(pdf):
    """width/height of a PDF's first page, via pdfinfo or PyMuPDF."""
    try:
        import fitz
        p = fitz.open(pdf)[0].rect
        return p.width / p.height
    except Exception:
        return 1.5


def main(path):
    s = open(path).read()
    body = s[s.index(r'\maketitle'):s.index(r'\begin{acknowledgments}')]
    figdir = os.path.join(os.path.dirname(path), 'figures')
    figs = 0.0; figlist = []
    for m in re.finditer(r'\\begin\{(figure\*?)\}(.*?)\\end\{\1\}', body, re.S):
        wide = m.group(1).endswith('*')
        g = re.search(r'includegraphics\[[^\]]*\]\{([^}]*)\}', m.group(2))
        name = g.group(1) if g else '?'
        pdf = None
        for d in (figdir, os.path.join(figdir, 'bkgfit')):
            for ext in ('.pdf', '.png'):
                if os.path.exists(os.path.join(d, name + ext)):
                    pdf = os.path.join(d, name + ext)
        ar = aspect(pdf) if pdf and pdf.endswith('.pdf') else 1.5
        w = (300 / (0.5 * ar) + 40) if wide else (150 / ar + 20)
        figs += w; figlist.append((name, wide, ar, w))
    tabs = 0.0; tablist = []
    for m in re.finditer(r'\\begin\{(table\*?)\}(.*?)\\end\{\1\}', body, re.S):
        wide = m.group(1).endswith('*')
        rows = m.group(2).count(r'\\')
        w = (26 + 13 * rows) if wide else (13 + 6.5 * rows)
        tabs += w; tablist.append((wide, rows, w))
    # strip floats, comments, commands; count what remains
    txt = re.sub(r'\\begin\{(figure\*?|table\*?)\}.*?\\end\{\1\}', ' ', body, flags=re.S)
    txt = re.sub(r'(?m)%.*$', '', txt)
    txt = re.sub(r'\\(cite|ref|label)\{[^}]*\}', ' ', txt)
    txt = re.sub(r'\\[a-zA-Z@]+\*?(\[[^\]]*\])?', ' ', txt)
    txt = re.sub(r'[{}$&_^~\\]', ' ', txt)
    words = len([w for w in txt.split() if re.search(r'[A-Za-z0-9]', w)])
    total = words + figs + tabs
    print(f'text words        {words}')
    for n, wide, ar, w in figlist:
        print(f'figure {n:32s} {"double" if wide else "single"}-col aspect {ar:.2f} -> {w:5.0f}')
    for wide, rows, w in tablist:
        print(f'table  {"double" if wide else "single"}-col {rows} lines -> {w:5.0f}')
    print(f'TOTAL             {total:.0f} / {LIMIT}   ({"OK" if total <= LIMIT else "OVER by %d" % (total - LIMIT)})')


if __name__ == '__main__':
    main(sys.argv[1])
