#!/usr/bin/env bash
set -euo pipefail
DATA_DIR="${1:-/path/to/raw/data}"
python src/verify_raw_inputs.py --repo . --data-dir "$DATA_DIR"
python src/phase2b_crispr.py --data-dir "$DATA_DIR" --out-dir results/phase2b
python src/phase3a_pharmacogenomics.py --data-dir "$DATA_DIR" --out-dir results/phase3a_pharmacogenomics
python src/audit_external_outputs.py --repo .
python src/audit_canonical.py --repo .
python src/make_figures.py
