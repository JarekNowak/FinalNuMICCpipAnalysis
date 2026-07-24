#!/usr/bin/env python3
# Minimal, note-specific Markdown -> LaTeX converter for NuMIInternalNoteCC1pi.md.
# Keeps Unicode (rendered by xelatex + DejaVu), converts headers/tables/lists/
# code/blockquotes, and embeds the five comparison figures in section 6.3.
import re, sys

SRC="NuMIInternalNoteCC1pi.md"
OUT="NuMIInternalNoteCC1pi.tex"
FIGDIR="/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer/unfold_output"
FIGS=[("plot_costhetamu_0.pdf", r"Unfolded d$\sigma$/d\,cos\,$\theta_\mu$ (fake data), four generators + tune + fake-data truth, with data/MC ratio panel."),
      ("plot_pmu_0.pdf",       r"Unfolded d$\sigma$/d$p_\mu$."),
      ("plot_ppi_0.pdf",       r"Unfolded d$\sigma$/d$p_\pi$."),
      ("plot_costhetapi_0.pdf",r"Unfolded d$\sigma$/d\,cos\,$\theta_\pi$."),
      ("plot_thetamupi_0.pdf", r"Unfolded d$\sigma$/d$\theta_{\mu\pi}$.")]

# Selection- and unfolding-stage diagnostics (cos theta_mu), embedded at end of §3.
SELFIGS=[("selection_costhmu_step1_reco_spectrum.pdf", r"Selected reco cos\,$\theta_\mu$ spectrum: data (GENIE fake data) vs MC signal + MC/EXT background."),
         ("selection_costhmu_step2_bkgd_subtraction.pdf", r"Background-subtracted reco spectrum with the MC signal overlaid (input to unfolding)."),
         ("selection_costhmu_step3_smearing_matrix.pdf", r"Smearceptance (response) matrix, reco vs true cos\,$\theta_\mu$."),
         ("selection_costhmu_step4_efficiency.pdf", r"Selection efficiency vs true cos\,$\theta_\mu$.")]

# Detector-variation overlays (reco top, true bottom), one per observable, at end of §4.
DVFIGS=[("detvar_pmu.pdf",     r"Detector variations, $p_\mu$: reco spectrum (top) and true signal (bottom), all nine samples POT-scaled to the CV."),
        ("detvar_ppi.pdf",     r"Detector variations, $p_\pi$."),
        ("detvar_costhmu.pdf", r"Detector variations, cos\,$\theta_\mu$."),
        ("detvar_costhpi.pdf", r"Detector variations, cos\,$\theta_\pi$."),
        ("detvar_thmupi.pdf",  r"Detector variations, $\theta_{\mu\pi}$.")]

def esc(s):
    # escape LaTeX specials in normal text (Unicode kept as-is)
    s=s.replace('\\',r'\textbackslash{}')
    for a,b in [('&',r'\&'),('%',r'\%'),('$',r'\$'),('#',r'\#'),
                ('_',r'\_'),('{',r'\{'),('}',r'\}'),('~',r'\textasciitilde{}'),
                ('^',r'\textasciicircum{}')]:
        s=s.replace(a,b)
    return s

def inline(s):
    # protect inline code first
    codes=[]
    def stash(m):
        codes.append(m.group(1)); return f"\x00{len(codes)-1}\x00"
    s=re.sub(r'`([^`]+)`',stash,s)
    s=esc(s)
    # bold / italic (after escaping; markers are ASCII * which esc leaves alone)
    s=re.sub(r'\*\*([^*]+)\*\*',r'\\textbf{\1}',s)
    s=re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)',r'\\emph{\1}',s)
    # restore code spans as \texttt with their own escaping
    def unstash(m):
        c=codes[int(m.group(1))]
        c=c.replace('\\',r'\textbackslash{}')
        for a,b in [('&',r'\&'),('%',r'\%'),('$',r'\$'),('#',r'\#'),('_',r'\_'),
                    ('{',r'\{'),('}',r'\}'),('~',r'\textasciitilde{}'),('^',r'\textasciicircum{}')]:
            c=c.replace(a,b)
        return r'\texttt{'+c+'}'
    s=re.sub('\x00(\d+)\x00',unstash,s)
    return s

def figures_block(figset=FIGS,width=0.86):
    out=[]
    for fn,cap in figset:
        out.append(r'\begin{figure}[htbp]\centering')
        out.append(rf'\includegraphics[width={width}\textwidth]{{{FIGDIR}/{fn}}}')
        out.append(rf'\caption{{{cap}}}')
        out.append(r'\end{figure}')
    return "\n".join(out)

lines=open(SRC,encoding='utf-8').read().split('\n')
body=[]; i=0; n=len(lines)
def flush_para(buf):
    if buf:
        body.append(inline(' '.join(buf))); body.append('')
para=[]
while i<n:
    ln=lines[i]
    # fenced code
    if ln.startswith('```'):
        flush_para(para); para=[]
        i+=1; code=[]
        while i<n and not lines[i].startswith('```'):
            code.append(lines[i]); i+=1
        i+=1
        body.append(r'\begin{verbatim}'); body.extend(code); body.append(r'\end{verbatim}'); body.append('')
        continue
    # headers
    m=re.match(r'^(#{1,3})\s+(.*)$',ln)
    if m:
        flush_para(para); para=[]
        lvl=len(m.group(1)); txt=inline(m.group(2))
        cmd={1:'section',2:'section',3:'subsection'}[lvl]
        # emit comparison figures right before section 6.4
        if m.group(2).strip().startswith('6.4'):
            body.append(figures_block()); body.append('')
        # emit selection/stage figures at the end of section 3 (before "4 ...")
        if re.match(r'^4\s',m.group(2).strip()):
            body.append(figures_block(SELFIGS,0.72)); body.append('')
        # emit detector-variation overlays at the end of section 4 (before "5 ...")
        if re.match(r'^5\s',m.group(2).strip()):
            body.append(figures_block(DVFIGS,0.80)); body.append('')
        body.append(rf'\{cmd}*{{{txt}}}'); body.append('')
        i+=1; continue
    # horizontal rule
    if re.match(r'^---+\s*$',ln):
        flush_para(para); para=[]
        body.append(r'\vspace{0.5em}\hrule\vspace{0.5em}'); body.append('')
        i+=1; continue
    # table (line with | and next line is separator)
    if ln.strip().startswith('|') and i+1<n and re.match(r'^\s*\|?[\s:|-]+\|',lines[i+1]):
        flush_para(para); para=[]
        header=[c.strip() for c in ln.strip().strip('|').split('|')]
        i+=2; rows=[]
        while i<n and lines[i].strip().startswith('|'):
            rows.append([c.strip() for c in lines[i].strip().strip('|').split('|')]); i+=1
        ncol=len(header)
        colspec='|'+'|'.join([r'>{\raggedright\arraybackslash}X']*ncol)+'|'
        body.append(r'\begin{tabularx}{\textwidth}{'+colspec+'}')
        body.append(r'\hline')
        body.append(' & '.join(r'\textbf{'+inline(h)+'}' for h in header)+r' \\')
        body.append(r'\hline')
        for r in rows:
            r=(r+['']*ncol)[:ncol]
            body.append(' & '.join(inline(c) for c in r)+r' \\')
        body.append(r'\hline')
        body.append(r'\end{tabularx}'); body.append('')
        continue
    # blockquote
    if ln.startswith('>'):
        flush_para(para); para=[]
        q=[]
        while i<n and lines[i].startswith('>'):
            q.append(lines[i].lstrip('>').strip()); i+=1
        # join into paragraphs on blank lines
        qtext=[]; cur=[]
        for qq in q:
            if qq=='':
                if cur: qtext.append(inline(' '.join(cur))); cur=[]
            else: cur.append(qq)
        if cur: qtext.append(inline(' '.join(cur)))
        body.append(r'\begin{quote}\itshape')
        body.append('\n\n'.join(qtext))
        body.append(r'\end{quote}'); body.append('')
        continue
    # lists (unordered - / ordered N.)
    mm=re.match(r'^(\s*)([-*]|\d+\.)\s+(.*)$',ln)
    if mm:
        flush_para(para); para=[]
        ordered=bool(re.match(r'\d+\.',mm.group(2)))
        env='enumerate' if ordered else 'itemize'
        body.append(rf'\begin{{{env}}}')
        while i<n:
            m2=re.match(r'^(\s*)([-*]|\d+\.)\s+(.*)$',lines[i])
            if not m2:
                # continuation line (indented, non-empty, not a new block)
                if lines[i].strip() and lines[i].startswith('  ') and not lines[i].strip().startswith('|'):
                    body[-1]=body[-1]+' '+inline(lines[i].strip()); i+=1; continue
                break
            body.append(r'\item '+inline(m2.group(3))); i+=1
        body.append(rf'\end{{{env}}}'); body.append('')
        continue
    # blank
    if ln.strip()=='':
        flush_para(para); para=[]
        i+=1; continue
    para.append(ln.strip()); i+=1
flush_para(para)

# title = first H1
title="NuMI CC1$\\pi$ Internal Note"
for l in lines:
    if l.startswith('# '):
        title=inline(l[2:].strip()); break

preamble=r'''\documentclass[11pt]{article}
\usepackage{fontspec}
\setmainfont{DejaVu Serif}[Scale=0.92]
\setsansfont{DejaVu Sans}[Scale=0.92]
\setmonofont{DejaVu Sans Mono}[Scale=0.85]
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{tabularx}
\usepackage{longtable}
\usepackage{parskip}
\usepackage[hidelinks]{hyperref}
\setlength{\emergencystretch}{3em}
\renewcommand{\arraystretch}{1.2}
\begin{document}
'''
# drop the first H1 from body (we render it as \title-ish heading) -> our loop already made it a \section*.
tex=preamble+'\n'.join(body)+'\n\\end{document}\n'
open(OUT,'w',encoding='utf-8').write(tex)
print("wrote",OUT,"(",len(body),"blocks )")
