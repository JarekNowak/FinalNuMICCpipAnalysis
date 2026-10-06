# report/

Documents, figures, tables and the data release of the MicroBooNE NuMI CC1π± cross-section analysis.

## Where things are

| Folder | Contents |
|---|---|
| `notes/` | The analysis documents, which cross-reference each other and are built together: `analysis_note` (inclusive), `proton_tagged_note`, `technical_supplement`, `change_log`, `approval_request`, and `cr_data_note` (control regions, kept separate from the three analysis documents). Their shared inputs: `shared/`, `glossary.tex`, `common_preamble.tex`, `status.tsv` (status of every result; read by `tools/inventory.py` and `tools/check_documents.py`) and `beam_on_held_out.{tex,json}` (passages held out until the control regions are reported). |
| `papers/` | PRD (inclusive) and PRL (proton-tagged) drafts, with `paper_bib.tex` and `paper_macros.tex` |
| `slides/` | `analysis_slides`, `review_panel_slides`, `collab_talk_imperial_2026` |
| `other_notes/` | `sigma0_note`, `multipion_note`, `newobs_note`, `unfolding_note` |
| `planning/` | Plans, reviews and the TODO list |
| `results/` | `current_results.tsv` and `current_results_1p_new.tsv` (written by `tools/make_results_tsv.py`, compared with the notes by `tools/check_tables.py`), `closure_summary.tsv` (`tools/closure_tables.C`), `normalisation_manifest.tsv` (`xsec_analyzer/norm_manifest.py`) and other result tables |
| `figures/`, `tables/` | Figures and generated LaTeX tables, used by every document |
| `data_release/` | Released cross sections, covariances, `A_C` matrices and model curves |
| `tools/` | Table, figure and check scripts: `python3 report/tools/<name>.py` from the repository root |
| `macros/` | ROOT figure macros (`root -l -b -q report/macros/<name>.C` from the repository root) and `eps2pdf.sh` |
| `releases/` | Packages sent to the Editorial Board (`release_v1.6/`, `release_v1.6.1/`) and the data-release zip of 2026-09-29 |
| `archive/` | Earlier Overleaf exports and superseded files |
| `references/`, `internalDocs/` | Cited papers, and internal notes of other analyses |

## Building the analysis documents

```bash
cd report/notes
for pass in 1 2; do
  for d in proton_tagged_note analysis_note technical_supplement change_log approval_request cr_data_note; do
    pdflatex -interaction=nonstopmode $d.tex > /dev/null
  done
done
grep -c '^!' *.log                     # errors: 0 expected
grep -c 'LaTeX Warning: Reference' *.log   # undefined references: 0 expected
```

The second pass resolves the cross-references between the documents (`xr`, prefixes `N-`, `P-`, `S-`, `CR-`, `AR-`).
Figures are read from `../figures/` and generated tables from `../tables/`. The papers, slides and other notes build
in their own folders in the same way; the papers also need `bibtex` between `pdflatex` runs.

## Checks after a change

```bash
python3 report/tools/check_tables.py       # note tables against results/current_results.tsv
python3 report/tools/check_documents.py    # status labels, prerequisites, terminology
python3 xsec_analyzer/check_staleness.py   # figures and tables newer than their sources
python3 xsec_analyzer/norm_manifest.py     # normalisation of every configuration
```
