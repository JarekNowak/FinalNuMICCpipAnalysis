"""check_documents.py -- cross-document consistency checks of the analysis package (review of 2026-09-27).

The analysis note, the proton-tagged note and the technical supplement must agree on the numbers that
define the measurement, use one status vocabulary, and take the shared statements from one source.
This refuses the package (exit 1) when
  1. a superseded exposure, flux or trigger value appears in prose (they drifted in earlier releases);
  2. the exposures, fluxes or target count appear with a value other than the adopted one;
  3. a status label other than the four of shared/status_vocab.tex is used ("proposed secondary", ...);
  4. a document does not \\input the shared decision box, prerequisites or validation ladder, or carries
     its own copy of them;
  5. a prerequisite identifier (U1..., P1..., V1, R1) is used that shared/prerequisites.tex does not define;
  6. the status column of data_release/index_extractions.tsv disagrees with report/status.tsv;
  7. an old term remains in prose (report/tools/terminology.py --check).
    python3 report/tools/check_documents.py
"""
import os, re, sys, csv, subprocess

HERE = os.path.dirname(os.path.abspath(__file__)); REP = os.path.dirname(HERE)
DOCS = ['analysis_note.tex', 'proton_tagged_note.tex', 'technical_supplement.tex']
SHARED = {'decision_box': ['analysis_note.tex', 'proton_tagged_note.tex', 'technical_supplement.tex'],
          'prerequisites': ['analysis_note.tex', 'proton_tagged_note.tex', 'technical_supplement.tex'],
          'validation_ladder': ['analysis_note.tex', 'proton_tagged_note.tex', 'technical_supplement.tex']}
STALE = [r'3\.283\s*\\times\s*10\^\{20\}', r'9\\?,?846\\?,?635', r'8\.857\s*\\times', r'6\.60865', r'19\.939']
ADOPTED = {  # quantity: (regex capturing the value, adopted value)
    'FHC exposure': (r'FHC[^.]{0,60}?\$([0-9.]+)\\times10\^\{20\}\$~POT', {'7.766'}),
    'RHC exposure': (r'RHC[^.]{0,60}?\$([0-9.]+)\\times10\^\{20\}\$~POT', {'11.082'}),
    'FHC flux': (r'FHC\s*\$?([0-9.]{7})\$?[,;]\s*RHC', {'6.81159'}),
    'target count': (r'N_\{?\\?mathrm\{Ar\}\}?\s*=\s*([0-9.]+)\\times10\^\{29\}', {'8.710', '8.71'}),
}
OLD_STATUS = [r'proposed secondary', r'secondary scope', r'validation product(?!s? \()', r'\bwithdrawn\b']


def strip(t):
    return re.sub(r'(?<!\\)%[^\n]*', '', t)


def main():
    bad = []
    txt = {d: strip(open(os.path.join(REP, d)).read()) for d in DOCS}
    for d, t in txt.items():
        for rx in STALE:
            for m in re.finditer(rx, t):
                near = t[max(0, m.start() - 400):m.end() + 400].lower()
                if 'open-trigger' in near or 'open trigger' in near:
                    continue          # quoted deliberately as the open-trigger figure it is
                bad.append(f'{d}: superseded value "{m.group(0)}"')
        for name, (rx, ok) in ADOPTED.items():
            for m in re.finditer(rx, t):
                if m.group(1) not in ok:
                    bad.append(f'{d}: {name} quoted as {m.group(1)} (adopted {sorted(ok)})')
        prose = re.sub(r'\\(?:label|[a-z]*ref)\{[^}]*\}(?:\{[^}]*\})?', ' ', t)   # label names are not wording
        for rx in OLD_STATUS:
            for m in re.finditer(rx, prose, flags=re.I):
                ctx = prose[max(0, m.start() - 40):m.end() + 20].replace('\n', ' ')
                bad.append(f'{d}: status wording outside the vocabulary: "...{ctx}..."')
    for f, users in SHARED.items():
        body = strip(open(os.path.join(REP, 'shared', f + '.tex')).read())
        key = re.search(r'\\label\{([^}]*)\}', body)
        for d in users:
            if f'\\input{{shared/{f}}}' not in txt[d]:
                bad.append(f'{d}: does not \\input shared/{f}.tex')
            if key and f'\\label{{{key.group(1)}}}' in txt[d]:
                bad.append(f'{d}: carries its own copy of {key.group(1)}')
    if '\\input{shared/decision_box}' not in txt['analysis_note.tex'].split('\\section{Measurement definition}')[0]:
        bad.append('analysis_note.tex: the decision box is not at the beginning')
    pre = open(os.path.join(REP, 'shared', 'prerequisites.tex')).read()
    defined = set(re.findall(r'^([UPVR]\d+) &', pre, flags=re.M))
    used = set()
    for d, t in txt.items():
        used |= set(re.findall(r'\b([UPVR]\d)\b(?=[\s,;.)\-]|--)', t))
    for l in open(os.path.join(REP, 'status.tsv')):
        if not l.startswith('#'):
            for rng in re.findall(r'([UPVR])(\d)--[UPVR]?(\d)', l):
                used |= {f'{rng[0]}{i}' for i in range(int(rng[1]), int(rng[2]) + 1)}
            used |= set(re.findall(r'\b([UPVR]\d)\b', l))
    for u in sorted(used - defined):
        bad.append(f'prerequisite identifier {u} is used but not defined in shared/prerequisites.tex')
    lab = {'primary': 'Primary result', 'secondary': 'Secondary result proposed for approval',
           'validation': 'Validation-only product', 'notreported': 'Not reported'}
    lines = [l for l in open(os.path.join(REP, 'status.tsv')) if l.strip() and not l.startswith('#')]
    st = {}
    for r in csv.DictReader(lines, delimiter='\t'):
        for k in r['observables'].split(','):
            if k not in ('-', 'total'):
                st[(r['family'], k)] = r
    for l in open(os.path.join(REP, 'data_release', 'index_extractions.tsv')):
        if l.startswith('#') or l.startswith('tag\t'):
            continue
        f = l.rstrip('\n').split('\t'); fam, cfg, obs, status = f[1], f[2], f[3], f[-1]
        r = st.get((fam, obs))
        if r is None:
            bad.append(f'index: {f[0]} has no row in status.tsv'); continue
        col = {'FHC5': 'FHC', 'RHCFULL': 'RHC', 'COMB': 'COMB'}[cfg]
        want = lab[r['status']] if r[col] == 'Y' or r['status'] == 'notreported' else lab['notreported']
        if status != want:
            bad.append(f'index: {f[0]} status "{status}", status.tsv gives "{want}"')
    t = subprocess.run([sys.executable, os.path.join(HERE, 'terminology.py'), '--check'] + [os.path.join(REP, d) for d in DOCS],
                       capture_output=True, text=True)
    if t.returncode:
        bad += ['terminology: ' + x for x in t.stdout.split('\n') if x]
    if bad:
        print(f'DOCUMENT CONSISTENCY FAILURES ({len(bad)}):'); [print('  ' + b) for b in bad]; sys.exit(1)
    print(f'documents consistent: {len(DOCS)} documents, {len(defined)} prerequisites, '
          f'{len(st)} released observable keys with a status')


if __name__ == '__main__':
    main()
