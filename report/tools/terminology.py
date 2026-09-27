"""terminology.py -- apply the standard terminology of the documents (review of 2026-09-27, item 6.4) to
prose, leaving labels, references, citations, file names, \\texttt content and graphics untouched.

  fake data / fake-data  -> pseudo-data         (formal prose; figure stamps keep "CV PSEUDO-DATA")
  sideband(s)            -> control region(s)   (in prose; labels such as sec:sidebands are unchanged)

    python3 report/tools/terminology.py file.tex [...]      rewrite in place and report the counts
    python3 report/tools/terminology.py --check file.tex    exit 1 if any old term remains in prose
"""
import re, sys

PROTECT = re.compile(r'(\\(?:label|ref|eqref|pageref|cite|input|include|includegraphics|texttt|url|href|'
                     r'noteref|suppref|ptref|crref|inclref|blindfig)\*?(?:\[[^\]]*\])?(?:\{[^{}]*\})+'
                     r'|%[^\n]*)')
RULES = [
    (re.compile(r'\bFake-data\b'), 'Pseudo-data'), (re.compile(r'\bfake-data\b'), 'pseudo-data'),
    (re.compile(r'\bFake(\s+)data\b'), 'Pseudo-data'), (re.compile(r'\bfake(\s+)data\b'), 'pseudo-data'),
    (re.compile(r'\bSidebands\b'), 'Control regions'), (re.compile(r'\bsidebands\b'), 'control regions'),
    (re.compile(r'\bSideband\b'), 'Control region'), (re.compile(r'\bsideband\b'), 'control region'),
]


def rewrite(text):
    parts = PROTECT.split(text); n = 0
    for i in range(0, len(parts), 2):          # even indices: unprotected prose
        for rx, rep in RULES:
            parts[i], k = rx.subn(rep, parts[i]); n += k
    return ''.join(parts), n


def main():
    args = sys.argv[1:]
    check = '--check' in args
    files = [a for a in args if a != '--check']
    bad = 0
    for f in files:
        t = open(f).read(); new, n = rewrite(t)
        if check:
            if n:
                print(f'{f}: {n} old terms in prose'); bad += n
        else:
            open(f, 'w').write(new); print(f'{f}: {n} replacements')
    sys.exit(1 if (check and bad) else 0)


if __name__ == '__main__':
    main()
