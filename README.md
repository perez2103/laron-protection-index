# LPI reproducibility repository

Reproducible computational workflows for the manuscript:

**A Laron inspired molecular state shows attenuated AKT mTOR signaling without targetable cancer vulnerability**

> **Interpretation note.** `Laron Protection Index (LPI)` is the historical project name. In the manuscript it is treated as a **Laron-inspired signaling-state score**, not a validated clinical cancer-protection biomarker.

## What this repository reproduces

1. **Phase 2B — DepMap CRISPR**
   - frozen 20-gene LPI
   - 914 solid-tumor cell lines / 20 lineages
   - 152-line BRCA/LUAD/COADREAD discovery set
   - 762-line non-overlapping validation set
   - 688-line lineage-held-out validation set
   - 17,787 Gene Effect features with >=90% pan-solid coverage
   - continuous models adjusted for lineage and proliferation
   - LPI stability, within-lineage permutations and expression-matched random signatures

2. **Phase 3A — Pharmacogenomics**
   - PRISM 24Q2
   - GDSC1 / GDSC2
   - same frozen LPI and proliferation covariate
   - predefined IGF1R / PI3K / AKT / mTOR family annotations
   - global and within-family FDR
   - cross-platform replication table

3. **Phase 3A — iLINCS**
   - exact UP/DOWN signature
   - genetic (LIB_6) and chemical (LIB_5) perturbational connectivity
   - leave-one-landmark-out sensitivity analysis

4. **Phase 3B — TCGA PanCanAtlas**
   - 9,632 primary solid tumors / 31 cancer types
   - proliferation endpoint
   - ABSOLUTE-purity sensitivity
   - RPPA AKT-mTOR phosphosignaling composite
   - random-effects meta-analysis
   - exploratory overall survival kept outside the primary gate

## Key audited conclusions

- **CRISPR final validation gate: 0 genes.**
- **Pharmacogenomic replication gate: negative.** PRISM has 0 global and 0 pathway-family FDR signals; GDSC1 contains one isolated within-family signal (AS605240) that does not survive global FDR and has no PRISM/GDSC2 replication; GDSC2 has 0 family FDR signals.
- **iLINCS perturbational connectivity: positive and node-specific.** GHR/AKT1 genetic perturbations and multiple PI3K/mTOR compounds are LPI-concordant; selected PIK3CB/PDPK1 perturbations are anti-concordant.
- **TCGA proliferation: negative pan-cancer result.**
- **TCGA AKT-mTOR RPPA: inverse association with LPI after purity adjustment.**
- **No docking/target nomination is supported.**

## Reproducibility audit correction

The audit discovered an important file-management issue in the historical working directory: the files named `LPI_Phase2B_discovery_continuous.csv` and `LPI_Phase2B_validation_nonoverlap.csv` were stale exports from an earlier parser run and contained **140** and **726** models, respectively. The raw data themselves were complete.

The analysis was re-run from the original DepMap 24Q2 inputs. The canonical repository outputs correctly contain **152 discovery** and **762 non-overlapping validation** models. The pan-solid output (914 lines) reproduced the historical full-run coefficients to numerical precision, and the scientific conclusion remains unchanged: discovery FDR contains zero genes and the final gate contains zero candidates.

The fixed-seed stochastic controls were also re-run with seed `20260915`. They give:

- within-lineage max-|t| empirical P = **0.056478**
- expression-matched random-signature empirical P = **0.960396**

These replace earlier exploratory-seed values while preserving the same negative conclusion.

See `audit/AUDIT_REPORT.md` and `audit/numerical_audit.csv`.

## Repository layout

```text
config/                         Frozen LPI and analysis definitions
data/                           Acquisition instructions + SHA-256 manifest
environment/                    Python dependency lock list
src/                            Local reproducible analysis scripts
notebooks/                      External-query notebooks (iLINCS / TCGA)
results/phase2b/                Canonical re-run Phase 2B tables
results/phase3a_pharmacogenomics/
results/phase3a_ilincs/         Processed API outputs used by the paper
results/phase3b_tcga/           Processed PanCanAtlas outputs used by the paper
audit/                          Numerical audit and claim traceability
figures/                        Figures regenerated from canonical outputs
manuscript/                     Audited manuscript/source-data snapshot
```

## Environment

Tested audit environment:

- Python 3.13.5
- NumPy 2.3.5
- pandas 2.2.3
- SciPy 1.17.0
- statsmodels 0.14.6
- Matplotlib 3.10.8
- openpyxl 3.1.5

Install:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
```

## Raw data

Do not commit the large third-party raw datasets to GitHub. Download the release-specific files listed in `data/README.md`, place them together in a local directory, and verify them against `data/raw_manifest.csv`.

## Run locally

```bash
# Verify exact raw files first
python src/verify_raw_inputs.py --repo . --data-dir /path/to/raw/data

# Phase 2B
python src/phase2b_crispr.py \
  --data-dir /path/to/raw/data \
  --out-dir results/phase2b

# Phase 3A pharmacogenomics (slower because the GDSC XLSX files are streamed)
python src/phase3a_pharmacogenomics.py \
  --data-dir /path/to/raw/data \
  --out-dir results/phase3a_pharmacogenomics

# Recompute checks on processed iLINCS and TCGA outputs
python src/audit_external_outputs.py --repo .

# Canonical release audit (does not require historical working files)
python src/audit_canonical.py --repo .

# Optional forensic comparison against historical working outputs, if available
python src/audit_outputs.py --repo . --legacy-dir /path/to/historical/working/files

# Rebuild figures
python src/make_figures.py
```

Or:

```bash
make phase2b DATA=/path/to/raw/data
make phase3a DATA=/path/to/raw/data
make audit
```

## External API workflows

`notebooks/01_iLINCS_API_v2.ipynb` and `notebooks/02_iLINCS_leave_one_landmark_out.ipynb` query iLINCS. `notebooks/03_TCGA_PanCanAtlas_v3.ipynb` retrieves the open PanCanAtlas resources and produces the TCGA result bundle.

For a paper-exact rerun, compare fresh API/download outputs to the versioned processed tables rather than silently replacing them.

## Audit standard

Deterministic local models are considered reproduced when coefficients, standard errors, test statistics, P values and FDR values agree within `1e-5` absolute tolerance. External-query layers are audited by recalculating the reported summaries from the versioned processed source tables.

## Before public GitHub/Zenodo release

- author metadata are complete for sole author Juan José Pérez Cervera (ORCID 0000-0002-2657-6695; no institutional affiliation; corresponding email `septiembre2103@hotmail.com`);
- MIT code license is included; third-party data remain under source-resource terms;
- create the GitHub repository and update `repository-code` in `CITATION.cff`;
- create the tagged release `v1.0.0`;
- archive the release in Zenodo and add the DOI to the manuscript Code availability section.

## Release

This repository snapshot is frozen as **v1.0.0**. See `RELEASE_NOTES_v1.0.0.md`, `PRE_RELEASE_CHECKLIST.md` and `GITHUB_ZENODO_RELEASE_GUIDE.md`. Author metadata are complete for sole author Juan José Pérez Cervera, including corresponding email, no institutional affiliation, no specific funding, no competing interests, manuscript approval and accountability. The public GitHub repository is `https://github.com/perez2103/laron-protection-index`. Only the Zenodo DOI remains to be inserted after archival.

## Integrity

`data/raw_manifest.csv` contains exact SHA-256 hashes for all raw/downloaded inputs used in the audit. `audit/repository_manifest_sha256.csv` contains hashes for the repository release itself.


## Generative-AI assistance

OpenAI ChatGPT (GPT-5.6 Sol) was used as a computational and writing assistant during code development/checking, reproducibility packaging, and manuscript drafting/revision. The human author reviewed the work and retains responsibility for its accuracy and integrity. See `AI_USE_DISCLOSURE.md`.


## Canonical manuscript files

The current author-approved submission drafts are:

- `manuscript/LPI_Communications_Biology_manuscript_AUTHOR_FINAL.docx`
- `manuscript/LPI_Communications_Biology_cover_letter_AUTHOR_FINAL.docx`
- `manuscript/LPI_Supplementary_Information_AUTHOR_FINAL.docx`

Final author: **Juan José Pérez Cervera** (ORCID `0000-0002-2657-6695`), sole author, no institutional affiliation, corresponding email `septiembre2103@hotmail.com`, no specific funding and no declared competing interests. The author approved the manuscript and accepts responsibility for its accuracy and integrity. The final GitHub URL is `https://github.com/perez2103/laron-protection-index`. Only the Zenodo DOI remains to be supplied before journal submission.
