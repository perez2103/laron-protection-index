# GitHub + Zenodo release guide

This repository is scientifically and administratively frozen as **v1.0.0**. Author, correspondence, funding, conflict-of-interest and accountability metadata are complete. The public GitHub URL is now `https://github.com/perez2103/laron-protection-index`. Only the Zenodo DOI remains to be generated/inserted.

## 1. Generate the two remaining public identifiers

Author, correspondence, funding, conflict-of-interest and accountability metadata are already complete. Do not alter them unless the authorship record genuinely changes.

After the public repository is created:

- update `repository-code` in `CITATION.cff` with the final GitHub URL;
- after Zenodo mints the release DOI, insert that DOI into:
  - `manuscript/LPI_Communications_Biology_manuscript_AUTHOR_FINAL.docx` (Code availability),
  - `README.md`, and
  - optionally the `identifiers` field of `CITATION.cff`.

Do **not** change `config/lpi_definition.json`, `config/analysis_plan.json`, canonical result tables, statistical thresholds, fixed seed, or audit files after tagging v1.0.0 without incrementing the version.

## 2. Create the GitHub repository

Repository:

`https://github.com/perez2103/laron-protection-index`

Large third-party raw datasets must **not** be committed. They are excluded by `.gitignore` and are documented in `data/README.md` and `data/raw_manifest.csv`.

## 3. Validate the public checkout

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

## 4. Tag the release

```bash
git tag -a v1.0.0 -m "Audited manuscript reproducibility release"
git push origin v1.0.0
```

Create a GitHub Release from tag `v1.0.0` and use the text in `RELEASE_NOTES_v1.0.0.md`.

## 5. Archive in Zenodo

Recommended route:

1. Sign in to Zenodo with GitHub.
2. Enable the GitHub repository in Zenodo's GitHub integration.
3. Create/publish the GitHub `v1.0.0` release.
4. Zenodo will archive the release and mint a version DOI.
5. Copy the DOI into:
   - `CITATION.cff` (optional `identifiers` field),
   - the manuscript Code availability statement,
   - the GitHub README.

If using a manual Zenodo upload instead, upload the exact frozen v1.0.0 archive and use `.zenodo.json` as the metadata template.

## 6. Freeze the manuscript-linked version

Once the DOI has been minted, do not overwrite v1.0.0. Corrections should be released as v1.0.1 or later. Cite the exact version DOI submitted with the manuscript.

## 7. What reviewers should be able to verify

Without proprietary/restricted inputs, reviewers can:

- inspect every frozen definition and validation gate;
- reproduce figures from canonical processed outputs;
- recalculate iLINCS and TCGA reported summaries from versioned source tables;
- inspect the exact API/download notebooks;
- verify all repository and raw-input hashes.

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

Still required after Zenodo archival: Zenodo DOI.
