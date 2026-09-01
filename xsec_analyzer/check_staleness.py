#!/usr/bin/env python3
"""check_staleness.py -- refuse to ship an artefact older than the extraction it depicts.

Three times in this analysis a downstream product silently lagged a fix: figures after
the central-value weighting fix, the A_C matrices after the same fix, and figures again
after the duplicated-A_C fix. Each time the tables and text were correct and the pictures
were not, and nothing complained -- the scripts reported success in aggregate.

This compares every released artefact and every result figure against the closure sidecar
of the extraction it comes from, and fails if the artefact is older.

  python3 check_staleness.py           # exit 1 if anything is stale
  python3 check_staleness.py --list    # also list what is current
"""
import os, sys, glob, re, time

REPO    = os.path.dirname(os.path.abspath(__file__))
REPORT  = os.path.abspath(os.path.join(REPO, "..", "report"))
RELEASE = os.path.join(REPORT, "data_release")
FIGS    = os.path.join(REPORT, "figures")
PROC    = "/data/uboone/processed"

TAG2SIDECAR = {}   # release tag -> sidecar path
def sidecar_for(tag):
    """incl_FHC5_pmu -> closure_hists_xsec_FHC5_pmu.root ; 1p_* -> ..._ccpi1p_*"""
    fam, rest = tag.split("_", 1)
    stem = f"ccpi1p_{rest}" if fam == "1p" else rest
    return os.path.join(PROC, f"closure_hists_xsec_{stem}.root")

stale, missing, ok = [], [], []

# Artefacts produced inside a single unfold do not share a timestamp: the covariance is
# dumped a second or two before the closure sidecar. Without a tolerance every covariance
# file reads as stale, which would make the check noise. Real staleness in this analysis
# has always been hours or days, so a few minutes separates the two cleanly.
SAME_RUN_TOLERANCE_S = 300

def check(path, src, what):
    if not os.path.exists(src):
        missing.append((path, src, what)); return
    if os.path.getmtime(path) < os.path.getmtime(src) - SAME_RUN_TOLERANCE_S:
        stale.append((path, src, what,
                      (os.path.getmtime(src)-os.path.getmtime(path))/3600.))
    else:
        ok.append((path, what))

# ---- released artefacts: each maps to exactly one extraction -----------------------
for p in glob.glob(os.path.join(RELEASE, "A_C_*.tsv")):
    check(p, sidecar_for(os.path.basename(p)[4:-4]), "A_C")
for p in glob.glob(os.path.join(RELEASE, "curves_*.tsv")):
    check(p, sidecar_for(os.path.basename(p)[7:-4]), "curves")
for d in glob.glob(os.path.join(RELEASE, "cov", "*")):
    if not os.path.isdir(d): continue
    src = sidecar_for(os.path.basename(d))
    for f in glob.glob(os.path.join(d, "*.txt")):
        check(f, src, "covariance")

# ---- result figures: map by the config/observable in the filename ------------------
# Only figures that actually depend on an extraction are checked; selection- and
# detector-level figures do not and are skipped rather than guessed at.
RESULT_PREFIX = ("dsigma", "ppi2bin", "systbreak", "wstep", "wsmear", "fw_xsec")
CFG = {"FHC5":"FHC5", "RHCFULL":"RHCFULL", "COMB":"COMB", "FHC":"FHC5", "RHC":"RHCFULL"}
unmapped = []
for p in sorted(glob.glob(os.path.join(FIGS, "*"))):
    b = os.path.basename(p)
    if not b.startswith(RESULT_PREFIX): continue
    cfg = next((v for k, v in CFG.items() if re.search(rf"[_.]{k}[_.]", b)), None)
    if not cfg:
        unmapped.append(b); continue
    fam = "1p" if "_1p" in b or b.startswith("wstep") or b.startswith("wsmear") else "incl"
    # the newest sidecar for that configuration is the right bound: a figure drawn from
    # any extraction of it must postdate the most recent one
    cands = glob.glob(os.path.join(PROC, f"closure_hists_xsec_{'ccpi1p_' if fam=='1p' else ''}{cfg}_*.root"))
    if not cands:
        unmapped.append(b); continue
    check(p, max(cands, key=os.path.getmtime), f"figure/{fam}")

print(f"  checked {len(ok)+len(stale)} artefacts")
print(f"    current : {len(ok)}")
print(f"    STALE   : {len(stale)}")
if missing: print(f"    no source found: {len(missing)}")
if unmapped: print(f"    figures not mapped to an extraction (skipped): {len(unmapped)}")

if stale:
    print("\n  STALE -- older than the extraction they depict:")
    for path, src, what, hrs in sorted(stale, key=lambda x: -x[3])[:25]:
        print(f"    {what:<12} {os.path.basename(path):<44} {hrs:6.1f} h older")
    if len(stale) > 25: print(f"    ... and {len(stale)-25} more")

if "--list" in sys.argv:
    print("\n  current:")
    for path, what in ok[:20]:
        print(f"    {what:<12} {os.path.basename(path)}")

sys.exit(1 if stale else 0)
