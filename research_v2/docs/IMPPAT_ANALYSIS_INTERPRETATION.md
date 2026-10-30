# NatDrug-AI IMPPAT screening: independent quality check

**Date:** 30 October 2026. **Scope:** analysis of the uploaded `NatDrugAI_IMPPAT_Repaired_Results.zip` and executed resume Jupyter notebook. All findings are from submitted output files; no external publications have been independently verified.

## Completion status

- IMPPAT preparation: completed. 18,314 chemical input records, 18,046 unique valid phytochemicals, 16 structures excluded due to CHEMBL203 training overlap.
- In-silico EGFR screening: completed for **18,030** structures using `RF_fingerprint_descriptors`. There are **5,246 scaffold groups including empty/acyclic scaffolds** according to the output manifest.
- Source-measured human EGFR: 57 exact IC50 measurement entries across 29 phytochemical IDs in audit; 13 structural nonoverlap candidates, but original experiments and articles have **not** been independently cross-checked.
- Repeated RF/MLP training from previous experiments was not part of this repaired ZIP; the notebook marks Experiments 1 and 4 as preserved.

## Screening distribution and applicability cautions

| Statistic | Value |
|---|---:|
| Screened molecules | 18,030 |
| Median predicted pIC50 | 4.9774 |
| Median maximum training Tanimoto | 0.2170 |
| Molecules with similarity < 0.4 | 16642 (92.30%) |
| Molecules with similarity >= 0.4 | 1388 (7.70%) |
| Missing/acyclic scaffold values in CSV | 2787 (15.46%) |
| Molecules flagged with scaffold represented in training | 4837 |
| Unique non-null flagged scaffold strings | 37 |

**Interpretation:** The manuscript should not report “4,837 scaffolds shared with training.” That field counts *molecules* flagged `Scaffold_in_training=True`, not unique scaffold strings. Missing scaffold strings appear as blank/NaN on CSV import and need an explicit empty-scaffold convention. The 0.4 threshold is not externally calibrated.

### Highest predicted activity (unrestricted; interpretation risky)

| Name | Predicted pIC50 | Predicted IC50 (nM) | Maximum Tanimoto |
|---|---:|---:|---:|
| Coleon A | 6.6032 | 249 | 0.165 |
| Changrolin | 6.4633 | 344 | 0.500 |
| Cinchonine | 6.1840 | 655 | 0.250 |
| Cinchonidine | 6.1840 | 655 | 0.250 |
| Quinidinone | 6.1757 | 667 | 0.290 |
| Quinidine | 6.1412 | 722 | 0.290 |
| Epiquinine | 6.1412 | 722 | 0.290 |
| Quinine | 6.1412 | 722 | 0.290 |
| Cinchonicinol | 6.1264 | 747 | 0.224 |
| Lobelanidine | 6.0838 | 825 | 0.352 |

### Higher-similarity screening subset (Tanimoto >= 0.4; still not validated)

| Name | Predicted pIC50 | Predicted IC50 (nM) | Maximum Tanimoto |
|---|---:|---:|---:|
| Changrolin | 6.4633 | 344 | 0.500 |
| Bavachalcone | 5.9019 | 1253 | 0.458 |
| Zeatin O-glucoside riboside | 5.8761 | 1330 | 0.516 |
| 7-Benzylaminooxazolo[5,4-d]pyrimidine | 5.8679 | 1356 | 0.548 |
| Broussochalcone A | 5.8661 | 1361 | 0.543 |
| 7-O-allylapigenin | 5.8581 | 1387 | 0.476 |
| 1-[2-(Benzyloxy)-6-hydroxyphenyl]ethan-1-one | 5.8564 | 1392 | 0.440 |
| Neobavachalcone | 5.8550 | 1396 | 0.458 |
| Kinetin riboside | 5.8501 | 1412 | 0.569 |
| Zeatin riboside | 5.8428 | 1436 | 0.589 |

## Thirteen provisional measured EGFR records

These cross-database candidate records **must not be claimed as independent external validation** until their primary source assay identifiers, papers, assay types and independence from ChEMBL training activities have been checked. The provided files have `Original_Reference = not stated` for all 13 candidates.

| Compound | Source | Measured IC50 (nM) | Measured pIC50 | Predicted pIC50 | Abs. error | Tanimoto |
|---|---|---:|---:|---:|---:|---:|
| Apigenin | BindingDB | 111000.0 | 3.955 | 5.534 | 1.579 | 0.500 |
| Kaempferol | BindingDB | 10400.0 | 4.983 | 5.718 | 0.735 | 0.784 |
| 4-Hydroxybenzoic acid | BindingDB | 72400.0 | 4.140 | 4.197 | 0.057 | 0.444 |
| biochanin A | BindingDB | 105000.0 | 3.979 | 5.726 | 1.747 | 0.750 |
| Acacetin | BindingDB | 140000.0 | 3.854 | 5.546 | 1.692 | 0.459 |
| Shikonin | NPASS | 139.2 | 6.856 | 4.934 | 1.922 | 0.222 |
| Prunetin | BindingDB | 14000.0 | 4.854 | 5.681 | 0.827 | 0.725 |
| Flavone | BindingDB | 224000.0 | 3.650 | 4.766 | 1.116 | 0.528 |
| Methyl 2,5-dihydroxycinnamate | BindingDB | 102000.0 | 3.991 | 4.871 | 0.880 | 0.405 |
| Sulforaphane | NPASS | 9000.0 | 5.046 | 4.576 | 0.470 | 0.136 |
| 1,5-Dihydroxyxanthone | NPASS | 900.0 | 6.046 | 4.807 | 1.239 | 0.361 |
| Orobol | BindingDB | 10400.0 | 4.983 | 5.821 | 0.838 | 0.784 |
| psi-Tectorigenin | BindingDB | 333.0 | 6.478 | 5.579 | 0.898 | 0.674 |

**Exploratory, non-validated diagnostic only (n=13)**: MAE = 1.077 pIC50, RMSE = 1.197 pIC50, R² = -0.388. The negative R² is a significant warning that the current model does not predict these particular measurements well. It is not a scientifically admissible external validation result until assay provenance is audited.

## Technical follow-ups before Molecular Informatics submission

1. Recover primary PubMed/DOI/BindingDB/NPASS assay identifiers for all 13 measured candidate records; compare assay source IDs, references, proteins, modalities, endpoints, and duplicate records against frozen ChEMBL203 training data. Confirm true external independence; eliminate duplicates or inconsistent endpoints.
2. Test model residuals versus similarity in **held-out scaffold folds** and report coverage–error results. Avoid treating 0.4 as an established reliability threshold.
3. Acknowledge that 92.3% of screened structures are under the unvalidated 0.4 similarity cutoff; do not propose ultra-low-similarity rank leaders (e.g., Coleon A) as confirmed EGFR inhibitors.
4. Verify whether the 4,837 “shared scaffold” number represents molecules, not scaffolds; standardize treatment and documentation of missing Murcko scaffolds. Fingerprint similarity can equal 1 even when the exact canonical structures differ.
5. Verify IMPPAT natural-product-origin and chemical identities for key-ranked molecules, and screen chemical alerts, false-positive risk, and physicochemical constraints before recommending experiments.
6. Report training-data curation and assay-sensitivity runs beside model benchmark outcomes. Never claim therapeutic efficacy, experimentally verified hits, or natural-product-specific model superiority from these predictions alone.

**Conclusion:** Experiment 3 in-silico screening completed; Experiment 2 has 13 *candidate* measured records but is not validated. The raw 13-record diagnostic suggests performance challenges that must be resolved or transparently disclosed before publication.