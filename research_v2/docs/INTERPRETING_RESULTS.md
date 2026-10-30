# How to interpret the numerical outputs

Results were read from the user's completed Jupyter analyses. They are not newly simulated data.

| Quantity | Measured status | Value | Notes |
|---|---|---:|---|
| CHEMBL203 curated structures | Source-backed | 3,619 | 4,991 raw IC50 records, assays pooled after curation |
| RF fingerprint+descriptor 5-fold scaffold R² | Internal CV | 0.6454 | MAE 0.6524 and RMSE 0.8580 pIC50 |
| ExtraTrees fingerprint+descriptor 5-fold scaffold R² | Internal CV | 0.6504 | Small numerical difference; not demonstrated significant |
| MLP 3-seed mean R² | Internal CV | about 0.5434 | Do not label external test |
| Natural products ranked | In silico only | 18,030 | IMPPAT 3.0; no measured inhibition implied |
| ChEMBL test empirical interval coverage | Held-out ChEMBL scaffolds | 86.5% | 90% nominal; may not generalize to natural products |
| Independently verified EGFR natural-product R² | Not available | — | Original assay/publication provenance missing |

**Scaffold interpretation:** a molecule without a ring-based Bemis–Murcko scaffold is `NO_MURCKO_FRAMEWORK`, not `SHARED_FRAMEWORK`. The old screening output mistakenly counted acyclic molecules as sharing a blank scaffold. The corrected code and `correct_scaffolds.py` do not alter pIC50 predictions.

**Model provenance:** model comparisons, robustness, and internal CV are based on available data; cross-validation performance does not guarantee prospective activity.
