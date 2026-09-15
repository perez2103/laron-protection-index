# LPI reproducibility repository

Reproducible computational workflows for the manuscript:

**A Laron inspired molecular state shows attenuated AKT mTOR signaling without targetable cancer vulnerability**

> **Interpretation note.** `Laron Protection Index (LPI)` is the historical project name. In the manuscript it is treated as a **Laron-inspired signaling-state score**, not a validated clinical cancer-protection biomarker.

## Main result

The study separates **molecular-state reproducibility** from **therapeutic vulnerability**:

- **CRISPR:** 0 genes pass the frozen discovery-to-validation gate.
- **PRISM/GDSC:** no replicated LPI-associated drug-sensitivity phenotype.
- **iLINCS:** selected GHR–PI3K–AKT–mTOR perturbations reproducibly induce an LPI-like transcriptomic state.
- **TCGA:** higher LPI is associated with lower purity-adjusted AKT–mTOR phosphosignaling, but not with a reproducible pan-cancer reduction in proliferation.
- **No docking or single-target nomination is supported.**

## Frozen LPI

Axis genes (14): `GHR, IGF1, IGF1R, IRS1, IRS2, PIK3CA, PIK3CB, PDPK1, AKT1, AKT2, MTOR, RPTOR, RPS6KB1, EIF4EBP1`.

FOXO-output genes (6): `BCL2L11, CDKN1B, GADD45A, SESN3, SOD2, CAT`.

`LPI = 0.5 × (-mean z(axis)) + 0.5 × (mean z(FOXO output))`, followed by standardization within lineage/cancer type.

The score was not trained on Laron syndrome patient samples and is not a clinical protection biomarker.

## Audited cohorts and gates

### Phase 2B — DepMap CRISPR

- DepMap Public 24Q2
- 914 solid-tumor cell lines / 20 lineages
- discovery: 152 BRCA/LUAD/COADREAD lines
- non-overlap validation: 762 lines
- lineage-held-out validation: 688 lines
- 17,787 Gene Effect features with >=90% pan-solid coverage
- continuous model: `GeneEffect ~ LPI + lineage + proliferation`
- final validation gate: **0 genes**

Fixed seed: `20260915`.

- within-lineage max-|t| empirical P = **0.056478**
- expression-matched random-signature empirical P = **0.960396**

### Phase 3A — pharmacogenomics

- PRISM: 637 LPI-scored models / 6,790 compounds; 0 global and 0 pathway-family FDR signals
- GDSC1: 448 models / 402 drugs; 0 global and 1 isolated pathway-family FDR signal
- GDSC2: 447 models / 295 drugs; 0 global and 0 pathway-family FDR signals

The isolated GDSC1 AS605240 association does not survive global FDR and lacks PRISM/GDSC2 replication.

### Phase 3A — iLINCS

The complete frozen query is six FOXO-output genes UP and 14 axis genes DOWN. Only 8/20 genes overlap the relevant L1000 landmark space, so leave-one-landmark-out analysis is a core robustness control.

Key robust connections include GHR and AKT1 genetic perturbations and multiple PI3K/mTOR compounds in the LPI-concordant direction, with CAT knockdown, TGX 221 and GSK-2334470 showing robust anti-concordance.

### Phase 3B — TCGA PanCanAtlas

- 9,632 primary solid tumors / 31 cancer types
- proliferation random-effects meta beta = **-0.012355**, P = **0.603393**
- purity-residualized proliferation beta = **-0.012819**, P = **0.589771**
- AKT–mTOR RPPA random-effects meta beta = **-0.038265**, P = **0.000688**, with 21/27 cancer types negative

## Audit correction

The reproducibility audit identified two historical working exports from an earlier incomplete-parser run (140 discovery and 726 validation models). They are **not canonical**. Re-reading the original DepMap 24Q2 CRISPR file with the versioned physical-line parser recovers the intended 152/762 cohorts. The scientific conclusion is unchanged: discovery FDR contains zero genes and the final validation gate contains zero candidates.

See [`audit/AUDIT_REPORT.md`](audit/AUDIT_REPORT.md).

## Public GitHub versus full archived snapshot

This GitHub repository is intentionally **lightweight and executable**. It contains:

- frozen definitions and analysis plan;
- analysis/audit scripts;
- exact raw-input SHA-256 manifest and acquisition instructions;
- iLINCS API notebooks and compact robustness source data;
- TCGA cancer-level source tables used to recalculate the manuscript meta-analyses;
- canonical result summaries and audit records.

Large third-party raw datasets are not redistributed. Large processed matrices, raw API payloads, binary manuscript files, rendered figures and the complete frozen result snapshot are retained in the **Zenodo v1.0.0 archive**. The DOI will be added here after archival.

## Repository layout

```text
config/                Frozen LPI and analysis definitions
data/                  Acquisition instructions + raw-input SHA-256 manifest
environment/           Python dependencies and audited environment
src/                   Analysis, audit and figure-generation code
notebooks/             External-query workflows
results/                Lightweight canonical summaries/source tables
audit/                  Audit report and generated audit tables
```

## Environment

Audited environment:

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

## Verify the lightweight checkout

The public source tables are sufficient to recalculate the key iLINCS robustness counts and the three TCGA meta-analyses:

```bash
python src/audit_external_outputs.py --repo .
python src/audit_canonical.py --repo .
python -m py_compile src/*.py
```

## Full raw-data rerun

Download the exact release-specific files described in [`data/README.md`](data/README.md), place them in a local directory, and verify them first:

```bash
python src/verify_raw_inputs.py --repo . --data-dir /path/to/raw/data
```

Then rerun the deterministic cell-line layers:

```bash
python src/phase2b_crispr.py \
  --data-dir /path/to/raw/data \
  --out-dir results/phase2b

python src/phase3a_pharmacogenomics.py \
  --data-dir /path/to/raw/data \
  --out-dir results/phase3a_pharmacogenomics
```

After those full tables are regenerated (or restored from the archived v1.0.0 snapshot), figures can be rebuilt with:

```bash
python src/make_figures.py
```

## External-query workflows

- `notebooks/01_iLINCS_API_v2.ipynb` — frozen full-signature LIB_6/LIB_5 query.
- `notebooks/02_iLINCS_leave_one_landmark_out.ipynb` — leave-one-landmark-out control.
- the exact TCGA PanCanAtlas v3 retrieval notebook is preserved in the complete v1.0.0 Zenodo archive; cancer-level paper source tables are versioned here for independent meta-analysis.

For a paper-exact rerun, compare fresh external-service outputs against the versioned release rather than silently replacing them.

## Raw-data integrity

`data/raw_manifest.csv` records SHA-256 hashes for every original raw/downloaded input used in the audit. These hashes identify the exact DepMap, PRISM, GDSC, iLINCS and TCGA inputs underlying v1.0.0.

## Release and citation

Version: **v1.0.0**  
Repository: https://github.com/perez2103/laron-protection-index

Release notes: [`RELEASE_NOTES_v1.0.0.md`](RELEASE_NOTES_v1.0.0.md)  
Release/deposition guide: [`GITHUB_ZENODO_RELEASE_GUIDE.md`](GITHUB_ZENODO_RELEASE_GUIDE.md)

The permanent **Zenodo DOI is pending**. Once minted, v1.0.0 will not be overwritten; subsequent changes will use a new semantic version.

## Author

**Juan José Pérez Cervera**  
ORCID: https://orcid.org/0000-0002-2657-6695  
Correspondence: `septiembre2103@hotmail.com`  
No institutional affiliation. No specific funding. No competing interests declared.

The author approved the manuscript and accepts responsibility for the accuracy, integrity and accountability of the work.

## Generative-AI assistance

OpenAI ChatGPT (GPT-5.6 Sol) was used as an interactive computational and writing assistant during code development/checking, reproducibility packaging and manuscript drafting/revision. The sole human author reviewed the analyses and outputs and retains responsibility for the work. See [`AI_USE_DISCLOSURE.md`](AI_USE_DISCLOSURE.md).

## License

Original repository code: **MIT License**. Third-party datasets remain subject to their source-resource terms and are not relicensed by this repository. See [`LICENSE_NOTICE.md`](LICENSE_NOTICE.md).
