# Reproducibility audit report

Audit date: **2026-09-15**  
Frozen seed for stochastic controls: **20260915**

## Executive result

**Overall numerical audit: PASS.** One file-management issue was identified and resolved before submission.

The original DepMap, PRISM and GDSC inputs were re-read from their downloaded release files. Deterministic Phase 2B and Phase 3A estimates reproduced the full-run outputs to numerical tolerance. iLINCS and TCGA summary statistics were independently recalculated from versioned processed source tables.

## Phase 2B correction identified during audit

Two historical working CSVs were stale exports from an earlier incomplete-parser run:

- discovery table: 140 models (stale)
- non-overlap validation table: 726 models (stale)

The raw CRISPR file was complete. The audited parser recovers the intended cohorts:

- pan-solid: **914** models / **20** lineages
- discovery: **152** models
- non-overlap validation: **762** models
- lineage-held-out validation: **688** models
- tested CRISPR genes: **17,787**

The correction does not change the conclusion: discovery contains **0 genes at FDR < 0.05** and the final validation gate contains **0 genes**.

The six pan-solid FDR associations remain SOX10, PSMB5, RASGRP3, TUBA1B, BLTP1 and CLINT1, but they are not promoted because the frozen discovery gate is empty.

NCOA4 in the audited pan-solid run: beta = **-0.0338481**, P = **9.627e-5**, FDR q = **0.111520**.

## Fixed-seed null controls

The stochastic controls were re-run with seed `20260915` and are now versioned in `src/phase2b_crispr.py`:

- observed maximum |t| = **4.714703**
- 300 within-lineage LPI permutations: empirical P = **0.056478**
- 100 expression-matched random 20-gene signatures: empirical P = **0.960396**

These replace exploratory-seed values while preserving the same negative inference.

## LPI stability

The full-cohort audit showed strong stability to gene omission and moderate outliers:

- leave-one-gene-out median Pearson r = **0.989459**
- leave-one-gene-out minimum r = **0.932897**
- winsorized expression r = **0.993031**
- within-lineage gene standardization r = **0.980930**
- median/MAD robust scaling r = **0.597700**

The last result defines an explicit normalization sensitivity.

## Phase 3A pharmacogenomics

Audited overlaps:

| Platform | Models | Global FDR <0.05 | Predefined-family FDR <0.05 |
|---|---:|---:|---:|
| PRISM | 637 | 0 | 0 |
| GDSC1 | 448 | 0 | 1 |
| GDSC2 | 447 | 0 | 0 |

The isolated GDSC1 family signal is AS605240: beta = **0.250278**, P = **0.000809**, family q = **0.033174**, global q = **0.196412**, n = **434**. It does not satisfy the frozen replication gate.

## Phase 3A iLINCS

The live iLINCS service is external, so the paper-exact query output is frozen in the archived release and the GitHub repository contains the query notebooks plus a compact key robustness table.

Audited leave-one-landmark-out examples:

- GHR: 9/9 positive; 7/9 FDR <0.05
- AKT1: 7/7 positive among evaluable conditions; 5/7 FDR <0.05
- CAT: 8/8 negative; 8/8 FDR <0.05
- KU0060648: 9/9 positive and significant
- LY-294002: 9/9 positive and significant
- TORIN-2: 8/8 positive and significant
- TGX 221: 9/9 negative and significant
- GSK-2334470: 5/5 negative and significant

This supports perturbational, node-specific state convergence; it is not evidence of drug sensitivity.

## Phase 3B TCGA

Random-effects meta-analysis was recalculated from cancer-specific source tables.

- primary solid tumors: **9,632**
- cancer types: **31**
- proliferation meta beta = **-0.012355**, P = **0.603393**, I2 = **76.05%**
- purity-residualized proliferation beta = **-0.012819**, P = **0.589771**
- AKT-mTOR RPPA meta beta = **-0.038265**, P = **0.000688**, I2 = **51.75%**
- random-signature proliferation empirical P = **0.574257**

Tumor purity remains an important limitation of the bulk-tissue LPI; the RPPA analysis is purity-adjusted.

## Scientific interpretation after audit

The audit supports the manuscript's constrained conclusion:

1. no reproducible single-gene dependency;
2. no replicated pharmacologic sensitivity phenotype;
3. robust perturbational transcriptomic convergence;
4. lower AKT-mTOR phosphosignaling in TCGA at higher LPI after purity adjustment;
5. no reproducible pan-cancer reduction in proliferation;
6. no target is eligible for docking.

## Public GitHub versus archived release

The public GitHub repository is intentionally lightweight: it contains code, frozen definitions, input hashes, notebooks/workflow documentation, compact verification tables and cancer-level source tables. Large third-party raw datasets are not redistributed. Large processed matrices, raw API payloads, binary manuscript files and the complete frozen result snapshot are retained in the versioned Zenodo archive for `v1.0.0`.
