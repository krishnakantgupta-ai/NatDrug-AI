# NatDrug-AI v2 — independent archive quality assurance (2026-10-30)

**Reviewed:** `NatDrugAI_v2_All_Outputs.zip`, 50 entries, ZIP integrity check passed.

## Verified numerical calculations

All eight original benchmark combinations in `legacy_benchmark/oof_predictions.csv` were recalculated using held-out predictions and reproduced the published CSV summary. For RF with fingerprints + descriptors (scaffold five-fold CV), N = 3,619, MAE = 0.652371, RMSE = 0.857995 and R² = 0.645379. For ExtraTrees (scaffold five-fold CV), N = 3,619, MAE = 0.649912, RMSE = 0.851860 and R² = 0.650432. All 3,619 observations occur once per model in the out-of-fold prediction files. The +0.00505 R² difference should not be claimed as a significant performance advantage without paired evaluation.

The separate, held-out ChEMBL scaffold set has N = 587, MAE = 0.741874, RMSE = 0.962321, R² = 0.597639. Exactly 508 of 587 labels (86.54%) fall inside the nominal 90% prediction intervals, calibrated with 378 ChEMBL observations. Interval half-width = 1.46334 pIC50. These intervals are **not calibrated for IMPPAT chemical space**.

The IMPPAT screening output has N = 18,030 unique canonical SMILES and 5,246 scaffold groups counting the empty/no-framework group. 93.78% of screened structures have maximum Tanimoto similarity below the exploratory cutoff of 0.4 to the smaller prediction-interval model training subset. This cutoff is not scientifically calibrated.

## Discovered scaffold annotation error

For the 18,030 screened molecules, 2,787 (15.46%) have an empty Bemis–Murcko scaffold, characteristic of acyclic structures. The original `Scaffold_seen_in_train` flag treated **all 2,787 as True** because the empty group was present in the model training set. This makes an apparent ring-framework overlap where no ring framework exists. The status must be **not applicable / no Murcko framework**, not shared.

The model-training-subset comparison yields:

| Framework status | Screened molecules |
| --- | ---: |
| Shared nonempty scaffold with the smaller reliability-model training subset | 563 |
| Nonempty scaffold not seen in that subset | 14,680 |
| No Murcko framework | 2,787 |

For comparison only, against the **full curated ChEMBL203 set** (not the smaller reliability-model training subset), 2,050 molecules have one of 37 shared *nonempty* scaffolds; 13,193 molecules have novel nonempty scaffolds; 2,787 have no ring framework. This full-set statistic must not be mislabeled as overlap with the reliability model's actual training subset.

`corrected_screening_scaffolds.csv` adds corrected scope-explicit status columns and preserves all original predictions and original flags for audit. The accompanying script can regenerate the correction from the original ZIP.

## Submission status

- **Completed:** internal ChEMBL scaffold and random CV, baseline comparison, three-seed MLP stability, curation sensitivity, model-comparison output, natural-product *in-silico* screening and held-out ChEMBL interval check.
- **Unresolved:** 13 IMPPAT candidate experimental EGFR IC50 rows have no primary paper DOI/PMID or original assay IDs. None is marked as verified independent. Accordingly `external_validation_metrics` is null. Do not call their provisional comparison external validation.
- **Reproducibility gap:** The 50-file final-results ZIP contains result tables, plots, and some software versions/checksums, but *not* the complete raw datasets and full pipeline scripts/model artifacts. The benchmark's `source_release` field still reads `TO BE RECORDED by investigator`. Release a reproducible source-and-data package and exact software environment before submission.
- **Chemical-space gap:** 15% of natural products have no Murcko ring framework and many of the others are far from training chemistry. Report this explicitly, and avoid describing uncalibrated intervals as reliable uncertainty for these molecules.
- **Model-development caution:** Do not choose/tune a final model using the 13 natural products until their independence has been established and a prospective evaluation plan is frozen.

**Overall:** Stronger exploratory software study with internally validated predictive benchmarks and a large-scale screening demonstration; **not yet a convincing external natural-product QSAR validation**.
