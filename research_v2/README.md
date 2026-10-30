# NatDrug-AI v2 — Research validation and natural-product screening

**Companion research package for the existing [NatDrug-AI repository](https://github.com/krishnakantgupta-ai/NatDrug-AI).**

This is an **additive folder**: it does not overwrite the original Streamlit application or claim to introduce a new molecular representation. It provides a source-aware EGFR benchmark, MLP stability, ExtraTrees comparison, applicability and uncertainty analysis, and a reproducible natural-product screening workflow.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r research_v2/code/requirements_v2.txt
python -m pip install notebook
jupyter notebook
```

Open `research_v2/notebooks/NatDrugAI_v2_Results_Viewer.ipynb` to inspect **included real outputs**. For the full study, first obtain the IMPPAT exports (see `sources/README.md`), then run `NatDrugAI_v2_Complete_Study_REQUIRES_IMPPAT_DATA.ipynb`.

## Study results (status on 2026-10-30)

| Evidence | Status |
|---|---|
| Measured CHEMBL203 RF, fingerprint+descriptors, five-fold scaffold CV: R² 0.6454 | Completed internal evaluation |
| ExtraTrees scaffold CV: R² 0.6504 | Completed internal evaluation; superiority not established |
| MLP three seeds: mean R² ~0.5434 | Completed internal evaluation |
| IMPPAT 18,030-structure in-silico virtual screen | Completed but **predictions only** |
| Held-out ChEMBL scaffold interval coverage 86.5% vs nominal 90% | Completed internal diagnostic |
| Independently verified natural-product EGFR validation | **Not established**; primary assay provenance needed |

## Contents

- `code/`: reproducible data-processing, benchmark, screening, validation-gate and report scripts, plus corrected Murcko scaffold handling.
- `data/`: ChEMBL IC50 source export; see CC BY-SA attribution and caveats.
- `notebooks/`: results viewer and full analysis (requires local IMPPAT files).
- `results/`: ChEMBL prediction files and model results; screening figures and aggregate manifests only (bulk IMPPAT screening records are omitted from public package).
- `docs/`: scientific quality audit, limitations, data attribution, reproducibility and upload checklist.
- `tests/`: scaffold status regression test.
- `sources/`: instructions for obtaining official bulk exports; these are **not automatically available from GitHub**.

## Scientific interpretation

Training molecules can overlap with externally collated assay reports even when their standardized structures are not identical. NPASS and BindingDB may quote the same underlying measurements. The 13 IMPPAT EGFR candidate values are **not a verified independent validation set**. Experimental follow-up is required before claiming biochemical inhibition.

The Tanimoto 0.4 cutoff is an *exploratory warning*, not a validated probability of correctness; predictions on structurally novel natural products should be treated cautiously.

Read `docs/INTERPRETING_RESULTS.md`, `docs/SCIENTIFIC_QA.md`, and `docs/DATA_SOURCES_AND_LICENSES.md` before submitting a manuscript.
