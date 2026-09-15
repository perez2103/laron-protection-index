# GitHub + Zenodo release guide

This repository is scientifically and administratively frozen as **v1.0.0**. Author, correspondence, funding, conflict-of-interest and accountability metadata are complete.

Public identifiers:

- GitHub: `https://github.com/perez2103/laron-protection-index`
- Zenodo DOI: `10.5281/zenodo.22771220`

## Current status

- GitHub repository: public and populated.
- Zenodo v1.0.0 archive: published under DOI `10.5281/zenodo.22771220`.
- Manuscript Code availability: updated to cite GitHub and Zenodo.
- Scientific definitions, thresholds, seed and canonical results remain frozen.

## Final GitHub release step

Create the GitHub tag/release `v1.0.0` from the manuscript-linked frozen state. Recommended release title:

`v1.0.0 — audited manuscript reproducibility release`

Use `RELEASE_NOTES_v1.0.0.md` as the release description.

Do **not** modify `config/lpi_definition.json`, `config/analysis_plan.json`, canonical result tables, statistical thresholds, fixed seed or audit outputs after tagging v1.0.0 without incrementing the version.

## Validate the public checkout

In a clean clone:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
pip install -r environment/requirements.txt
python src/audit_external_outputs.py --repo .
python src/audit_canonical.py --repo .
python -m py_compile src/*.py
```

For full local raw-data re-runs, place the release-specific third-party inputs in a local data directory and follow `README.md`.

## Archive policy

The manuscript-linked reproducibility archive is:

`https://doi.org/10.5281/zenodo.22771220`

Do not overwrite the scientific meaning of v1.0.0. Corrections after release should use v1.0.1 or later and, where appropriate, a new Zenodo version record.

## What reviewers can verify

Without redistributing third-party raw data, reviewers can:

- inspect every frozen definition and validation gate;
- reproduce figures from canonical processed outputs;
- recalculate iLINCS and TCGA reported summaries from versioned source tables;
- inspect API/download notebooks;
- verify repository and raw-input hashes.

With the publicly downloadable raw DepMap/PRISM/GDSC inputs, reviewers can additionally rerun the deterministic cell-line analyses from scratch.

## Final author metadata

- Sole author: Juan José Pérez Cervera
- ORCID: 0000-0002-2657-6695
- Corresponding email: septiembre2103@hotmail.com
- No institutional affiliation
- No specific funding
- No competing interests declared
- Manuscript approved 15 September 2026
- Author accepts responsibility for the accuracy and integrity of the work
