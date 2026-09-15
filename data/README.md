# Data acquisition and integrity

Raw third-party datasets are **not committed to this repository**. Place them in a local directory and pass that directory with `--data-dir` (or `DATA=/path/to/raw/data` in the Makefile).

The audit used the exact files and SHA-256 hashes listed in `raw_manifest.csv`. A rerun should first verify the hashes. If a source has reissued a file, use the archived/release-specific version rather than silently substituting a current release.

## DepMap / PRISM

Release: **DepMap Public 24Q2** for model metadata, expression and CRISPR Gene Effect, plus **Repurposing Public 24Q2** for PRISM.

Expected filenames:

- `Model(2).csv` (source name `Model.csv`)
- `OmicsExpressionProteinCodingGenesTPMLogp1(2).csv` (source name `OmicsExpressionProteinCodingGenesTPMLogp1.csv`)
- `CRISPRGeneEffect (1).csv` (source name `CRISPRGeneEffect.csv`)
- `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv`
- `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv`

Release archive: https://plus.figshare.com/articles/dataset/DepMap_24Q2_Public/25880521
DepMap portal: https://depmap.org/portal/data_page/

### Important parser note

The CRISPR file used in the project contains 1,150 physical data rows and 18,443 Gene Effect columns. One quoting anomaly can confuse a standards-compliant CSV parser and make many physical rows appear merged. The reproducibility code therefore reads the file by physical line, treats the first field as `ModelID`, and parses the remaining numeric comma-delimited fields directly. No values are imputed or repaired.

Both copies that were uploaded during the analysis (`CRISPRGeneEffect (1).csv` and `CRISPRGeneEffect(2).csv`) were bitwise identical in the original working session; the manifest records the canonical `(1)` copy.

## GDSC

Expected files:

- `GDSC1_fitted_dose_response_27Oct23.xlsx`
- `GDSC2_fitted_dose_response_27Oct23.xlsx`

Current access is via Cell Model Passports / GDSC resources:
https://cellmodelpassports.sanger.ac.uk/downloads

The scripts map `SANGER_MODEL_ID` to the DepMap `SangerModelID` field and use `LN_IC50` as the primary response, with `AUC` as a sensitivity endpoint.

## iLINCS

The repository contains the exact Colab notebooks used to query the public iLINCS API and the processed API outputs used in the paper. The API workflow follows the documented `SignatureMeta/upload`, `findConcordances` and `signatureEnrichment` routes.

Documentation/example API workflow:
https://github.com/uc-bd2k/ilincsAPI/blob/master/usingIlincsApis.Rmd

Because iLINCS is a live service, fresh API results can in principle change if its libraries are updated. For exact paper reproduction, use the processed outputs and hashes recorded in `raw_manifest.csv`; use the notebooks for a fresh validation query.

## TCGA PanCanAtlas

The TCGA notebook downloads open-access PanCanAtlas resources through GDC, including:

- batch-adjusted RNA-seq expression
- `TCGA-RPPA-pancan-clean.txt`
- ABSOLUTE purity/ploidy calls
- TCGA Clinical Data Resource

PanCanAtlas page:
https://gdc.cancer.gov/about-data/publications/pancanatlas

For exact paper reproduction without relying on future GDC changes, the processed Phase 3B outputs are versioned under `results/phase3b_tcga/` and audited from cancer-specific source tables.
