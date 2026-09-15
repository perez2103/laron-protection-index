PYTHON ?= python
DATA ?= /path/to/raw/data

verify:
	$(PYTHON) src/verify_raw_inputs.py --repo . --data-dir "$(DATA)"

phase2b:
	$(PYTHON) src/phase2b_crispr.py --data-dir "$(DATA)" --out-dir results/phase2b

phase3a:
	$(PYTHON) src/phase3a_pharmacogenomics.py --data-dir "$(DATA)" --out-dir results/phase3a_pharmacogenomics

audit-external:
	$(PYTHON) src/audit_external_outputs.py --repo .

audit:
	$(PYTHON) src/audit_external_outputs.py --repo .
	$(PYTHON) src/audit_canonical.py --repo .

audit-legacy:
	$(PYTHON) src/audit_outputs.py --repo . --legacy-dir "$(LEGACY)"

figures:
	$(PYTHON) src/make_figures.py

all: verify phase2b phase3a audit figures
