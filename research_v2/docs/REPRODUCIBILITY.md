# Reproducing the study

Prerequisites: Python 3.10+ with RDKit, NumPy, pandas, scikit-learn, scipy, matplotlib, Jupyter and optional PyTorch; see `code/requirements_v2.txt`. Tested previously across macOS and the report-generation patch, but **the GitHub upload package has not rerun every model from scratch**.

From a fresh clone of the existing NatDrug-AI repository, run:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r research_v2/code/requirements_v2.txt
python -m pip install notebook
jupyter notebook
```

Open `research_v2/notebooks/NatDrugAI_v2_Results_Viewer.ipynb` to review the included real outputs (works without the IMPPAT bulk tables).

For a full fresh rerun:
1. Download official IMPPAT chemical and bioactivity exports; place them as described in `sources/README.md`.
2. Open `research_v2/notebooks/NatDrugAI_v2_Complete_Study_REQUIRES_IMPPAT_DATA.ipynb`.
3. On Step 2, change the `FORCE_REBUILD_*` options if you need to recompute saved results. Rebuilding can be computationally intensive.
4. Generate `natural_products_scaffolds_corrected.csv` using `code/correct_scaffolds.py` if inspecting an older cached screening CSV.
5. Read the provenance gate status before interpreting any natural-product IC50 results.

For local screening source files, the separate local-only archive can be extracted into `research_v2/` to restore cached inputs, avoiding repeated expensive simulations or curation.

Optional automated test setup: copy `docs/GITHUB_ACTIONS_TEMPLATE.yml` to the repository root at `.github/workflows/natdrug_v2_tests.yml`. GitHub will only execute workflow files in the repository-root `.github/workflows/`.
