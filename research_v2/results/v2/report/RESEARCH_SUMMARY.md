# NatDrug-AI v2 — Evidence-first computational study report

Report generated from supplied measured datasets and machine-generated outputs. **Not a manuscript-ready external validation claim.**

## Source and reproducibility

- ChEMBL source file SHA256: `08eea3ccde0c01953f0a7bb46f863b256df480571f61edad4e64da6cc71ce042`
- IMPPAT bioactivity table SHA256: `6a8aeaa8d7350c2513decf52733da50fbe1917e88cd6f187299a8b88b007e33e`
- Human EGFR: ChEMBL target CHEMBL203, UniProt P00533.
- IMPPAT 3.0 release: 2026-09-30.

## Experiment 1: measured ChEMBL assay benchmark

| Split    | Model                      |   Compounds |   Scaffolds |   MAE_pIC50 |   RMSE_pIC50 |      R2 |
|:---------|:---------------------------|------------:|------------:|------------:|-------------:|--------:|
| random   | RF_fingerprint             |        3619 |        1204 |      0.5921 |       0.7762 |  0.7098 |
| random   | RF_fingerprint_descriptors |        3619 |        1204 |      0.5416 |       0.7347 |  0.7399 |
| random   | median                     |        3619 |        1204 |      1.1963 |       1.4476 | -0.0095 |
| random   | ridge                      |        3619 |        1204 |      0.963  |       1.4403 |  0.0007 |
| scaffold | RF_fingerprint             |        3619 |        1204 |      0.6859 |       0.8825 |  0.6249 |
| scaffold | RF_fingerprint_descriptors |        3619 |        1204 |      0.6524 |       0.858  |  0.6454 |
| scaffold | median                     |        3619 |        1204 |      1.2062 |       1.4611 | -0.0284 |
| scaffold | ridge                      |        3619 |        1204 |      1.3958 |       2.0717 | -1.0676 |

Numbers above are cached actual previous-Mac-run results, NOT newly computed by this report builder.

## Experiment 4: MLP stability

|   Seed |   Scaffolds |   Compounds |   MAE_pIC50 |   RMSE_pIC50 |     R2 |
|-------:|------------:|------------:|------------:|-------------:|-------:|
|     42 |        1204 |        3619 |      0.7482 |       0.9668 | 0.5497 |
|     52 |        1204 |        3619 |      0.7754 |       0.9954 | 0.5227 |
|     62 |        1204 |        3619 |      0.7302 |       0.9581 | 0.5578 |

## Alternative models

| Model              | Split    |    N |   MAE_pIC50 |   RMSE_pIC50 |     R2 |
|:-------------------|:---------|-----:|------------:|-------------:|-------:|
| ExtraTrees_fp_desc | scaffold | 3619 |      0.6499 |       0.8519 | 0.6504 |

## Experiment 3: natural-product virtual screening + uncertainty

- In-silico-screened IMPPAT structures: 18,030
- ChEMBL scaffold-held-out interval coverage (NOT natural-product external coverage): 0.865
- Interval half-width: 1.463 pIC50
- Screening structures below heuristic Tanimoto 0.4: 93.8%
- Intervals do not have guaranteed coverage on novel natural-product chemistry.

## Experiment 2: externally measured natural-product EGFR IC50

- Status: **BLOCKED_UNVERIFIED_SOURCE_PUBLICATION_ASSAY_METADATA**
- Candidate structures: 13
- Reason: IMPPAT Source names are not primary assay/publication IDs; ChEMBL training dataset has only assay IDs and no primary document IDs.

## Interpretation and submission blockers

1. Exact SMILES non-overlap does not rule out shared assay or publication sources.
2. IMPPAT NPASS / BindingDB pairs may contain the same reported experimental measurements.
3. Independent EGFR natural-product performance remains **unverified**; do not call 13-row diagnostic metrics an external validation.
4. Prediction intervals calibrated on ChEMBL scaffolds may under-cover structurally novel natural products.
5. Freeze all datasets, seed/split assignments, source release and environment before submission.
6. A negative natural-product diagnostic result must be transparently reported; model improvements cannot be chosen using an untouched external test.
