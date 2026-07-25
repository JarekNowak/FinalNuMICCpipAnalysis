#!/usr/bin/env python3
# Note-specific Markdown -> LaTeX converter for NuMIInternalNoteCC1pi.md.
# Keeps Unicode (rendered by xelatex + DejaVu); converts headers/tables/lists/
# code/blockquotes. Figures are inserted via explicit markers in the markdown:
#   <!--FIGSET:name-->     emits the named figure block (see FIGSETS below).
# Missing figure files are silently skipped so a partial plot set still compiles.
import re, os

SRC="NuMIInternalNoteCC1pi.md"
OUT="NuMIInternalNoteCC1pi.tex"
FIGDIR="/home/t2k/nowak/MicroBooNE/working_xsec_analyzer/xsec_analyzer/unfold_output"

OBS   = ["pmu","ppi","costhmu","costhpi","thmupi"]
OLAB  = {"pmu":r"$p_\mu$","ppi":r"$p_\pi$","costhmu":r"cos\,$\theta_\mu$",
         "costhpi":r"cos\,$\theta_\pi$","thmupi":r"$\theta_{\mu\pi}$"}
# comparison-plot basenames differ from the internal obs keys
CMP   = {"pmu":"pmu","ppi":"ppi","costhmu":"costhetamu","costhpi":"costhetapi","thmupi":"thetamupi"}

def sel(step, tmpl):   # per-observable selection-stage figure list
    return [(f"selection_{o}_{step}.pdf", tmpl.format(o=OLAB[o])) for o in OBS]

COVFIGS=[("XsecUnits_mat_table_cov_total.pdf",       r"Total covariance (systematic $+$ statistical) on the unfolded cos\,$\theta_\mu$ cross section."),
         ("XsecUnits_mat_table_cov_flux_total.pdf",  r"Flux (PPFX) covariance."),
         ("XsecUnits_mat_table_cov_xsec_total.pdf",  r"GENIE cross-section (UBGenie $+$ unisim) covariance."),
         ("XsecUnits_mat_table_cov_detVar_total.pdf",r"Detector-variation covariance."),
         ("XsecUnits_mat_table_cov_reint.pdf",       r"Hadron-reinteraction covariance."),
         ("XsecUnits_mat_table_cov_add_smear.pdf",   r"Additional smearing matrix $A_C$ (the Wiener-SVD extra smearing applied to predictions).")]

FIGSETS = {
 "comparison": ([ (f"plot_{CMP[o]}_0.pdf",
        (r"Unfolded d$\sigma$/d"+OLAB[o]+r" (fake data) on the \emph{optimised} binning (\S3.8): four generators (GENIE, GiBUU, NuWro, NEUT) $+$ MicroBooNE tune $+$ fake-data truth (all $A_C$-smeared), with data/MC ratio panel."
         if o=="pmu" else r"Unfolded d$\sigma$/d"+OLAB[o]+r" (optimised binning, four generators).")) for o in OBS], 0.84),
 "selstage":   (sel("step1_reco_spectrum",  r"Selected reco {o} spectrum: data (GENIE fake data) vs MC signal $+$ MC/EXT background.")[2:3]
              + sel("step2_bkgd_subtraction",r"Background-subtracted {o} reco spectrum with MC signal overlaid (unfolding input).")[2:3]
              + sel("step3_smearing_matrix", r"Smearceptance (response) matrix, reco vs true {o}.")[2:3]
              + sel("step4_efficiency",      r"Selection efficiency vs true {o}.")[2:3], 0.70),
 "reco_all":   (sel("step1_reco_spectrum",  r"Reco {o}: data (fake) vs MC signal $+$ background at the final cut."), 0.62),
 "resp_all":   (sel("step3_smearing_matrix",r"Response (smearceptance) matrix, reco vs true {o}."), 0.60),
 "eff_all":    (sel("step4_efficiency",     r"Selection efficiency vs true {o}."), 0.60),
 "detvar":     ([ (f"detvar_{o}.pdf",
        (r"Detector variations, "+OLAB[o]+r": reco spectrum (top) and true signal (bottom), all nine samples POT-scaled to the CV."
         if o=="pmu" else r"Detector variations, "+OLAB[o]+r".")) for o in OBS], 0.80),
 "cov":        (COVFIGS, 0.60),
 "systbreak":  ([ (f"syst_breakdown_{o}.pdf",
        r"Fractional systematic uncertainty by source, "+OLAB[o]+r": total (black), cross section (red), flux/beamline (blue), detector (green), reinteraction (magenta), MC$+$data stats (orange), POT$+$targets (grey).") for o in OBS], 0.62),
 "regmatrix":  ([("plot_regularization_matrix.pdf", r"Second-derivative regularisation matrix $C$ used in the Wiener-SVD (tridiagonal, Neumann boundary).")], 0.55),
 "mares":      ([("mares_band_costhmu.pdf", r"M$_A^{\mathrm{RES}}\pm1\sigma$ band on the GENIE CC1$\pi$ cos\,$\theta_\mu$ cross section (nominal black; $+1\sigma$ red, $-1\sigma$ blue). Asymmetric and Q$^2$-dependent (largest at backward angles).")], 0.68),
}

def esc(s):
    s=s.replace('\\',r'\textbackslash{}')
    for a,b in [('&',r'\&'),('%',r'\%'),('$',r'\$'),('#',r'\#'),
                ('_',r'\_'),('{',r'\{'),('}',r'\}'),('~',r'\textasciitilde{}'),
                ('^',r'\textasciicircum{}')]:
        s=s.replace(a,b)
    return s

def inline(s):
    codes=[]
    def stash(m):
        codes.append(m.group(1)); return f"\x00{len(codes)-1}\x00"
    s=re.sub(r'`([^`]+)`',stash,s)
    s=esc(s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'\\textbf{\1}',s)
    s=re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)',r'\\emph{\1}',s)
    def unstash(m):
        c=codes[int(m.group(1))]
        c=c.replace('\\',r'\textbackslash{}')
        for a,b in [('&',r'\&'),('%',r'\%'),('$',r'\$'),('#',r'\#'),('_',r'\_'),
                    ('{',r'\{'),('}',r'\}'),('~',r'\textasciitilde{}'),('^',r'\textasciicircum{}')]:
            c=c.replace(a,b)
        return r'\texttt{'+c+'}'
    s=re.sub('\x00(\d+)\x00',unstash,s)
    return s

def figures_block(name):
    figset,width=FIGSETS[name]
    out=[]
    for fn,cap in figset:
        if not os.path.exists(os.path.join(FIGDIR,fn)):
            continue                         # skip missing plots gracefully
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
    # figure marker
    mfig=re.match(r'^\s*<!--\s*FIGSET:(\w+)\s*-->\s*$',ln)
    if mfig:
        flush_para(para); para=[]
        name=mfig.group(1)
        if name in FIGSETS:
            body.append(figures_block(name)); body.append('')
        i+=1; continue
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
    m=re.match(r'^(#{1,4})\s+(.*)$',ln)
    if m:
        flush_para(para); para=[]
        lvl=len(m.group(1)); txt=inline(m.group(2))
        cmd={1:'section',2:'section',3:'subsection',4:'subsubsection'}[lvl]
        body.append(rf'\{cmd}*{{{txt}}}'); body.append('')
        i+=1; continue
    # horizontal rule
    if re.match(r'^---+\s*$',ln):
        flush_para(para); para=[]
        body.append(r'\vspace{0.5em}\hrule\vspace{0.5em}'); body.append('')
        i+=1; continue
    # table
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
    # lists
    mm=re.match(r'^(\s*)([-*]|\d+\.)\s+(.*)$',ln)
    if mm:
        flush_para(para); para=[]
        ordered=bool(re.match(r'\d+\.',mm.group(2)))
        env='enumerate' if ordered else 'itemize'
        body.append(rf'\begin{{{env}}}')
        while i<n:
            m2=re.match(r'^(\s*)([-*]|\d+\.)\s+(.*)$',lines[i])
            if not m2:
                if lines[i].strip() and lines[i].startswith('  ') and not lines[i].strip().startswith('|'):
                    body[-1]=body[-1]+' '+inline(lines[i].strip()); i+=1; continue
                break
            body.append(r'\item '+inline(m2.group(3))); i+=1
        body.append(rf'\end{{{env}}}'); body.append('')
        continue
    if ln.strip()=='':
        flush_para(para); para=[]
        i+=1; continue
    para.append(ln.strip()); i+=1
flush_para(para)

title="NuMI CC1$\\pi$ Internal Note"
for l in lines:
    if l.startswith('# '):
        title=inline(l[2:].strip()); break

preamble=r'''\documentclass[11pt]{article}
\usepackage{fontspec}
\setmainfont{DejaVu Serif}[Scale=0.92]
\setsansfont{DejaVu Sans}[Scale=0.92]
\setmonofont{DejaVu Sans Mono}[Scale=0.80]
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{tabularx}
\usepackage{longtable}
\usepackage{parskip}
\usepackage{float}
\usepackage[hidelinks]{hyperref}
\setlength{\emergencystretch}{3em}
\renewcommand{\arraystretch}{1.2}
\begin{document}
'''
tex=preamble+'\n'.join(body)+'\n\\end{document}\n'
open(OUT,'w',encoding='utf-8').write(tex)
print("wrote",OUT,"(",len(body),"blocks )")
