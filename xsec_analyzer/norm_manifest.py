#!/usr/bin/env python3
"""norm_manifest.py -- canonical normalisation audit for the CC1mu1pi NuMI analysis.

Every extraction divides by the same four numbers: the data POT, the integrated flux,
the argon target count, and the per-run MC exposure. Eight separate normalisation bugs
in this analysis's history came from one of those four drifting between configurations
without anything noticing -- an obsolete flux, an FHC flux reused for RHC, a run-total
POT duplicated across split files, a detVar exposure double-counted.

This script reads the authoritative sources (the xsec configs and the file_properties
tables -- not the note, which is downstream) and asserts that they agree, so a drift
becomes a failed check instead of a wrong cross section.

Usage:  python3 norm_manifest.py            # audit, write the TSV, exit 1 on failure
"""
import glob, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CFG  = os.path.join(HERE, "configs")
OUT  = os.path.abspath(os.path.join(HERE, "..", "report", "normalisation_manifest.tsv"))

# Argon nuclei in the fiducial volume, dead-wire slab excluded (analysis note Sec. 5).
N_AR = 8.710e29

# Reference values the note quotes. The manifest's job is to prove the configs still
# produce these, so a silent drift in either direction is caught.
EXPECT_POT  = {"fhc5": 8.857e20, "rhcfull": 1.1082e21, "comb": 1.9939e21}
EXPECT_FLUX = {"fhc5": 6.81159e-10, "rhcfull": 6.44646e-10, "comb": 6.60865e-10}

CONFIGS = ["fhc5", "rhcfull", "comb"]
failures, notes, rows = [], [], []

def fail(msg):
    failures.append(msg)

def note(msg):
    notes.append(msg)

def parse_fp(path):
    """file_properties: <file> <run> <type> [triggers] [pot]"""
    per_run_pot, per_run_files, types = defaultdict(float), defaultdict(list), defaultdict(int)
    if not os.path.exists(path):
        return None
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 3:
            continue
        fpath, run, ftype = f[0], f[1], f[2]
        types[ftype] += 1
        if ftype == "onBNB" and len(f) >= 5:
            # A run split across several files repeats its run-total POT in every file
            # (the processed/ convention), so take the DISTINCT value per run rather than
            # summing lines. Every current config lists one onBNB line per run, which makes
            # this a no-op -- but the two POT conventions in this analysis (run-total
            # duplicated in processed/, true per-file in processed/w/) are exactly what
            # produced an earlier normalisation bug, so the caller re-checks the
            # one-line-per-run assumption instead of quietly relying on it.
            per_run_pot[run] = max(per_run_pot[run], float(f[4]))
            per_run_files[run].append(fpath)
    return per_run_pot, per_run_files, types

for cfg in CONFIGS:
    # --- every observable's xsec config for this configuration must agree on flux ---
    xcs = sorted(glob.glob(os.path.join(CFG, f"ccpi*_xsec_config_numi_*_{cfg}.txt")))
    fluxes, fpfiles, unfolds = {}, {}, {}
    for xc in xcs:
        txt = open(xc).read()
        m = re.search(r"^\s*Flux\s+([0-9.eE+-]+)", txt, re.M)
        if m:
            fluxes[os.path.basename(xc)] = float(m.group(1))
        m = re.search(r"^\s*FPFile\s+(\S+)", txt, re.M)
        if m:
            fpfiles[os.path.basename(xc)] = m.group(1)
        m = re.search(r"^\s*Unfold\s+(.+)$", txt, re.M)
        if m:
            unfolds[os.path.basename(xc)] = m.group(1).strip()

    if not xcs:
        fail(f"[{cfg}] no xsec configs found")
        continue

    # CHECK 1: one flux per configuration, matching the authoritative value.
    #
    # Compared with a relative tolerance, not for exact string equality: some configs
    # spell the same number to a different number of digits (6.44646e-10 vs 6.446460e-10
    # are identical; comb has 6.60865e-10 vs 6.608653e-10, a 4.5e-7 relative difference).
    # Those are transcription artefacts worth knowing about but far below any level that
    # moves a cross section, and rewriting the configs to unify them would break the
    # byte-match between each config and the result it produced. So: tolerance for the
    # pass/fail verdict, and a separate informational note about the spellings.
    FLUX_RTOL = 1e-5
    uniq_flux = sorted(set(fluxes.values()))
    flux = uniq_flux[0] if uniq_flux else float("nan")
    if uniq_flux:
        spread = (max(uniq_flux) - min(uniq_flux)) / min(uniq_flux)
        if spread > FLUX_RTOL:
            fail(f"[{cfg}] FLUX DRIFT across observables (rel. spread {spread:.2e}): {uniq_flux}")
            for k, v in sorted(fluxes.items()):
                fail(f"    {k}: {v:.7g}")
        elif len(uniq_flux) > 1:
            note(f"[{cfg}] flux written {len(uniq_flux)} ways, max rel. difference "
                 f"{spread:.1e} (numerically equivalent): "
                 + ", ".join(f"{v:.7g}" for v in uniq_flux))
    if uniq_flux and abs(flux - EXPECT_FLUX[cfg]) / EXPECT_FLUX[cfg] > 1e-6:
        fail(f"[{cfg}] flux {flux:.6g} != expected {EXPECT_FLUX[cfg]:.6g}")

    # CHECK 2: one regularisation prescription per configuration (the review's "freeze").
    uniq_unf = sorted(set(unfolds.values()))
    if len(uniq_unf) != 1:
        fail(f"[{cfg}] UNFOLD PRESCRIPTION DRIFT: {uniq_unf}")

    # CHECK 3: POT, summed over DISTINCT per-run values.
    uniq_fp = sorted(set(fpfiles.values()))
    pot_by_fp = {}
    for fp in uniq_fp:
        p = fp if os.path.isabs(fp) else os.path.join(HERE, fp)
        parsed = parse_fp(p)
        if parsed is None:
            fail(f"[{cfg}] missing file_properties: {fp}")
            continue
        per_run_pot, per_run_files, types = parsed
        pot_by_fp[fp] = sum(per_run_pot.values())
        # every run needs BOTH onBNB and extBNB or the framework's map::at throws
        runs_on = set(per_run_pot)
        runs_ext = set()
        for line in open(p):
            f = line.split()
            if len(f) >= 3 and f[2] == "extBNB":
                runs_ext.add(f[1])
        # If a run ever spans several onBNB lines, the max-per-run rule above stops being
        # equivalent to the intended total and the convention must be settled explicitly.
        for r, fl in sorted(per_run_files.items()):
            if len(fl) > 1:
                fail(f"[{cfg}] {os.path.basename(fp)}: run {r} has {len(fl)} onBNB lines -- "
                     f"per-run POT aggregation is ambiguous, settle the convention")
        missing = runs_on - runs_ext
        if missing:
            fail(f"[{cfg}] {os.path.basename(fp)}: runs with onBNB but no extBNB: {sorted(missing)}")

    for fp, tot in pot_by_fp.items():
        if abs(tot - EXPECT_POT[cfg]) / EXPECT_POT[cfg] > 1e-3:
            fail(f"[{cfg}] {os.path.basename(fp)}: POT {tot:.5g} != expected {EXPECT_POT[cfg]:.5g}")

    pot = list(pot_by_fp.values())[0] if pot_by_fp else float("nan")
    rows.append((cfg, len(xcs), f"{flux:.6g}", f"{pot:.5g}", f"{N_AR:.4g}",
                 uniq_unf[0] if uniq_unf else "?", ";".join(os.path.basename(x) for x in uniq_fp)))

# ---------------------------------------------------------------------------------
# Systematics/file-list coupling.
#
# The proton-tagged detector systematic was absent from all 33 extractions because the
# detVar SAMPLES existed and were processed, but neither the file_properties table nor
# the systematics config referenced them. Neither file was wrong on its own; the two
# simply did not agree, and nothing compared them. A declared source with no files
# contributes nothing, and files with no declaration are silently ignored -- both are
# failures, and the second is the one that actually happened.
SYSTCALC = {"fhc5":  "ccpi_systcalc_numi.conf",
            "rhcfull":"ccpi_systcalc_numi.conf",
            "comb":  "ccpi_systcalc_numi.conf"}
SYSTCALC_1P = "ccpi1p_systcalc_numi.conf"

def declared_dv(conf):
    """DV sources declared in a systcalc config: '<name> DV <file_type>'."""
    out = set()
    p = os.path.join(CFG, conf)
    if not os.path.exists(p):
        return None
    for line in open(p):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) >= 3 and f[1] == "DV":
            out.add(f[2])
    return out

def available_dv(fp):
    """detVar file types present in a file_properties table."""
    out = {}
    p = fp if os.path.isabs(fp) else os.path.join(HERE, fp)
    if not os.path.exists(p):
        return None
    for line in open(p):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        f = line.split()
        if len(f) >= 3 and f[2].startswith("detVar"):
            out[f[2]] = out.get(f[2], 0) + 1
    return out

for cfg in CONFIGS:
    for conf, fpname, label in [
        (SYSTCALC[cfg],  f"file_properties_numi_{cfg}.txt",   f"{cfg} inclusive"),
        (SYSTCALC_1P,    f"file_properties_numi_{cfg}_w.txt", f"{cfg} proton-tagged")]:
        dec = declared_dv(conf)
        av  = available_dv(os.path.join(CFG, fpname))
        if dec is None or av is None:
            fail(f"[{label}] cannot read {conf} or {fpname}")
            continue
        # a declared source with no files contributes nothing
        for d in sorted(dec - set(av) - {"detVarCV"}):
            fail(f"[{label}] systematic '{d}' is declared but has NO files")
        # files with no declaration are silently ignored -- the failure that occurred
        for a in sorted(set(av) - dec - {"detVarCV"}):
            fail(f"[{label}] detVar sample '{a}' is present ({av[a]} file(s)) but NOT "
                 f"declared in {conf} -- it will be silently ignored")
        if dec and av and "detVarCV" not in av:
            fail(f"[{label}] detVar sources declared but no detVarCV reference sample")
        if dec and av:
            note(f"[{label}] {len(dec)} DV sources declared, all backed by files")

with open(OUT, "w") as fh:
    fh.write("config\tn_xsec_configs\tflux_cm2_per_POT\tdata_POT\tN_Ar\tunfold\tfile_properties\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

print(f"  wrote {OUT}")
print()
print(f"  {'config':<10} {'flux':<13} {'data POT':<12} {'unfold':<26} {'obs':>4}")
for r in rows:
    print(f"  {r[0]:<10} {r[2]:<13} {r[3]:<12} {r[5]:<26} {r[1]:>4}")
print()
if notes:
    print("  informational:")
    for n in notes:
        print(f"    {n}")
    print()
if failures:
    print(f"  FAILED {len(failures)} check(s):")
    for f in failures:
        print(f"    {f}")
    sys.exit(1)
print("  all normalisation checks PASSED")
